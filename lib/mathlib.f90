module mathlib
use lib_params, only:dp,pi

public :: invmati
public :: invreal33
public :: O3_to_axis_angle
public :: Dmatrix

contains
subroutine invmati(M,DIM)
      
    INTEGER           DET,M(3,3),TMP(3,3)
    INTEGER           DIM(3,3),DTMP(3,3),SUM,INA
    DOUBLE PRECISION  TES
!
      DET= M(1,1)*M(2,2)*M(3,3) + M(1,2)*M(2,3)*M(3,1)   &
          +M(1,3)*M(2,1)*M(3,2) - M(1,1)*M(2,3)*M(3,2)        &
          -M(1,2)*M(2,1)*M(3,3) - M(1,3)*M(2,2)*M(3,1)
!
      IF(ABS(DET).LT.1.E-4) STOP 'ERROR in invimat, DET=0'
!
      TMP(1,1)= (M(2,2)*M(3,3)-M(2,3)*M(3,2))
      TMP(1,2)=-(M(2,1)*M(3,3)-M(2,3)*M(3,1))
      TMP(1,3)= (M(2,1)*M(3,2)-M(2,2)*M(3,1))
!
      TMP(2,1)=-(M(1,2)*M(3,3)-M(1,3)*M(3,2))
      TMP(2,2)= (M(1,1)*M(3,3)-M(1,3)*M(3,1))
      TMP(2,3)=-(M(1,1)*M(3,2)-M(1,2)*M(3,1))
!
      TMP(3,1)= (M(1,2)*M(2,3)-M(1,3)*M(2,2))
      TMP(3,2)=-(M(1,1)*M(2,3)-M(1,3)*M(2,1))
      TMP(3,3)= (M(1,1)*M(2,2)-M(1,2)*M(2,1))
!
      DO 100 I=1,3
      DO 100 J=1,3
        INA=NINT(DBLE(TMP(J,I))/DBLE(DET))
        TES=DABS(DBLE(TMP(J,I))/DBLE(DET) - DBLE(INA))
        IF(TES.GT.1E-4) THEN
           STOP 'not integer inverse matrices'
        ELSE 
           DIM(I,J)=INA
        ENDIF
 100    CONTINUE

      SUM=0.D0
      DO 200  I=1,3
      DO 200  J=1,3
      DTMP(I,J)=0.D0
      DO 220 IJ=1,3
 220    DTMP(I,J)=DTMP(I,J)+DIM(I,IJ)*M(IJ,J)
      IF(I.EQ.J) SUM=SUM+ABS(DTMP(I,J)-1.0)
      IF(I.NE.J) SUM=SUM+ABS(DTMP(I,J)-0.0)
 200  CONTINUE
!     
      IF(SUM.GT.1.0E-6) STOP 'ERROR in INVMAT'
      RETURN

end subroutine invmati

subroutine  invreal33(mat,invmat)
implicit none
real(8) ,intent(in):: Mat(3,3)
real(8) ,intent(out):: Invmat(3,3)

real(8) :: pa(3),pb(3),pc(3),volumn
real(8) :: pra(3),prb(3),prc(3),t
  pa(:)=Mat(:,1); pb(:)=Mat(:,2); pc(:)=Mat(:,3)

  volumn =  pa(1)*(pb(2)*pc(3)-pb(3)*pc(2)) &
          + pa(2)*(pb(3)*pc(1)-pb(1)*pc(3)) &
          + pa(3)*(pb(1)*pc(2)-pb(2)*pc(1))

  t=1._8/volumn

  pra(1)=t*(pb(2)*pc(3)-pb(3)*pc(2))
  pra(2)=t*(pb(3)*pc(1)-pb(1)*pc(3))
  pra(3)=t*(pb(1)*pc(2)-pb(2)*pc(1))

  prb(1)=t*(pc(2)*pa(3)-pc(3)*pa(2))
  prb(2)=t*(pc(3)*pa(1)-pc(1)*pa(3))
  prb(3)=t*(pc(1)*pa(2)-pc(2)*pa(1))

  prc(1)=t*(pa(2)*pb(3)-pa(3)*pb(2))
  prc(2)=t*(pa(3)*pb(1)-pa(1)*pb(3))
  prc(3)=t*(pa(1)*pb(2)-pa(2)*pb(1))

 invMat(1,:) = pra(:); invMat(2,:) = prb(:); invMat(3,:) = prc(:)
end subroutine  invreal33

subroutine Dmatrix(ih, ik, il, rnx, rny, rnz, rtheta, rphi, degree, twoja1, Dmat)
!!ref    : https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
!!        s, px, py, pz, xy, yz, zx, x2-y2, 3z2-r2
!!        0,  1, -1,  0; -2, -1,   1,   2,  0
!!        fxyz, f5x3-xr2, f5y3-yr2, f5z3-zr2, fx(y2-z2), fy(z2-x2), fz(x2-y2)
!!TABLE I: http://journals.aps.org/prb/pdf/10.1103/PhysRevB.79.045107
  implicit none
  integer, intent(in), optional :: ih, ik, il
  real(dp), intent(in), optional :: rnx, rny, rnz
  real(dp), intent(in), optional :: rtheta, rphi
  real(dp), intent(in) ::  degree
  integer, intent(in) :: twoja1  ! equal to 2*J + 1
  ! Actually, this rotational Dmat should be a real matrix.
  complex(dp),dimension(twoja1,twoja1),intent(out) :: Dmat
  complex(dp),dimension(twoja1,twoja1)             :: Dtmp
  ! the key parameter for Dmatrix
  real(dp) :: nx, ny, nz, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4, nx5,ny5,nz5,nx6,ny6,nz6
  real(dp) ::  omega
  ! complex(dp),dimension(aint(J)*2+1,aint(J)*2+1) :: Jx,Jy,Jz
  ! integer :: m
  complex(dp),parameter:: ci = (0.d0,1.d0)
  real*8     ,parameter:: PI = 3.141592653589793238462643383279d0
  real(dp) :: tmp
  real(dp) :: sqrt3, sqrt15 
  real(dp) :: cosOmega_2, sinOmega_2
  real(dp) :: cos1omega, sin1omega, cos2omega, sin2omega, cos3omega, sin3omega, cos4omega, sin4omega
  real(dp) :: cos5omega, sin5omega, cos6omega, sin6omega
  real(dp) :: cos3omega_2, sin3omega_2, cos5omega_2, sin5omega_2 
  complex(dp) :: exp2omega, exp_2omega, exp4omega, exp_4omega, expomega_2, exp_omega_2
  complex(dp) :: exp3omega, exp_3omega, exp6omega, exp_6omega, expomega, exp_omega, exp5omega, exp_5omega  
  sqrt3 = dsqrt(3.0d0)
  sqrt15 = dsqrt(15.0d0)
  if(present(ih).and.present(ik).and.present(il))then
    tmp = dsqrt(dble(ih)*dble(ih) + dble(ik)*dble(ik) + dble(il)*dble(il))
    nx = dble(ih)/tmp
    ny = dble(ik)/tmp
    nz = dble(il)/tmp
  else if(present(rtheta).and.present(rphi))then
    nx = dsin(rtheta)*dcos(rphi)
    ny = dsin(rtheta)*dsin(rphi)
    nz = dcos(rtheta)
  else if(present(rnx).and.present(rny).and.present(rnz))then
    tmp = dsqrt(rnx*rnx + rny*rny + rnz*rnz )
    if(abs(tmp-1.d0).gt.1.d-6) then
      write(0,'(A,F16.10)')"WARNING!!! The nx, ny, nz input is not normailized: ",tmp
      write(0,*)"The program will normalize them now."
    end if
      nx = rnx/tmp
      ny = rny/tmp
      nz = rnz/tmp
  else
    write(*,*)" Please choose one of the three type of input: "
    write(*,*)" 1. h k l:"
    write(*,*)" 2. nx ny nz:"
    write(*,*)" 3. theta phi:"
  end if
  omega = dble(nint(degree))*PI/180.d0
  IF(twoja1 == 5) omega=omega/2.d0
  cosOmega_2 = dcos(omega/2.d0)
  sinOmega_2 = dsin(omega/2.d0)
  cos3omega_2 = dcos(3.d0*omega/2.d0)
  sin3omega_2 = dsin(3.d0*omega/2.d0)
  cos5omega_2 = dcos(5.d0*omega/2.d0)
  sin5omega_2 = dsin(5.d0*omega/2.d0)
  cos1omega = dcos(omega)
  sin1omega = dsin(omega)
  cos2omega = dcos(2.d0*omega)
  sin2omega = dsin(2.d0*omega)
  cos3omega = dcos(3.d0*omega)
  sin3omega = dsin(3.d0*omega)
  cos4omega = dcos(4.d0*omega)
  sin4omega = dsin(4.d0*omega)
  cos5omega = dcos(5.d0*omega)
  sin5omega = dsin(5.d0*omega)
  cos6omega = dcos(6.d0*omega)
  sin6omega = dsin(6.d0*omega)
  expomega = cmplx(cos1omega, sin1omega)
  exp_omega = cmplx(cos1omega, -sin1omega)
  exp2omega = cmplx(cos2omega,sin2omega)
  exp_2omega = cmplx(cos2omega,-sin2omega)
  exp4omega = cmplx(cos4omega,sin4omega)
  exp_4omega = cmplx(cos4omega,-sin4omega)
  exp3omega = cmplx(cos3omega,sin3omega)
  exp_3omega = cmplx(cos3omega,-sin3omega)
  exp5omega = cmplx(cos5omega, sin5omega)
  exp_5omega = cmplx(cos5omega, -sin5omega)
  exp6omega = cmplx(cos6omega,sin6omega)
  exp_6omega = cmplx(cos6omega,-sin6omega)
  expomega_2 = cmplx(cosOmega_2, sinOmega_2)
  exp_omega_2 = cmplx(cosOmega_2, -sinOmega_2)
  nx2 = nx**2
  ny2 = ny**2
  nz2 = nz**2
  nx3 = nx**3
  ny3 = ny**3
  nz3 = nz**3
  nx4 = nx**4
  ny4 = ny**4
  nz4 = nz**4
  nx5 = nx**5
  ny5 = ny**5
  nz5 = nz**5
  nx6 = nx**6
  ny6 = ny**6
  nz6 = nz**6
  if(twoja1 == 2) then
    Dmat(1,1) = cmplx(cosOmega_2, -nz*sinOmega_2)
    Dmat(2,2) = cmplx(cosOmega_2, nz*sinOmega_2)
    Dmat(1,2) = cmplx(-ny, -nx)*sinOmega_2
    Dmat(2,1) = cmplx(ny, -nx)*sinOmega_2
  else if(twoja1 == 3) then
    Dmat(1,1) = nx**2+(1-nx**2)*cos1omega
    Dmat(2,2) = ny**2+(1-ny**2)*cos1omega
    Dmat(3,3) = nz**2+(1-nz**2)*cos1omega
    Dmat(1,2) = nx*ny*(1-cos1omega) - nz*sin1omega
    Dmat(2,1) = nx*ny*(1-cos1omega) + nz*sin1omega
    Dmat(1,3) = nx*nz*(1-cos1omega) + ny*sin1omega
    Dmat(3,1) = nx*nz*(1-cos1omega) - ny*sin1omega
    Dmat(2,3) = ny*nz*(1-cos1omega) - nx*sin1omega
    Dmat(3,2) = ny*nz*(1-cos1omega) + nx*sin1omega
  else if(twoja1 == 5) then
    Dmat(1,1) = 3.d0*ny**2*nx**2 + &
                0.5d0*(exp_4omega + exp4omega)*(1-ny**2)*(1-nx**2) + &
                0.5d0*(exp_2omega + exp2omega)*(nz**2*(1-nz**2)+(nx**2-ny**2)**2)
    Dmat(2,2) = 3.d0*ny**2*nz**2 + &
                0.5d0*(exp_4omega + exp4omega)*(1-ny**2)*(1-nz**2) + &
                0.5d0*(exp_2omega + exp2omega)*(nx**2*(1-nx**2)+(ny**2-nz**2)**2)
    Dmat(3,3) = 3.d0*nz**2*nx**2 + &
                0.5d0*(exp_4omega + exp4omega)*(1-nz**2)*(1-nx**2) + &
                0.5d0*(exp_2omega + exp2omega)*(ny**2*(1-ny**2)+(nx**2-nz**2)**2)
    Dmat(4,4) = 0.75d0*(nx**2-ny**2)**2 + &
                0.125d0*(exp_4omega + exp4omega)*(4*nz**2+(nx**2-ny**2)**2) + &
                0.5d0*(exp_2omega + exp2omega)*(nz**2*(1-nz**2) + 4*ny**2*nx**2)
    Dmat(5,5) = 0.25d0*(1-3*nz**2)**2 + &
                1.5d0*(exp_2omega + exp2omega)*nz**2*(1-nz**2) + &
                0.375d0*(exp_4omega + exp4omega)*(1-nz**2)**2

    Dmat(1,2) = -2.d0*ci*sin1omega*(ci*ny**3*cos1omega + ci*ny*(1-ny**2)*cos3omega + &
                                      nx*nz*( 3.d0*ci*ny**2*sin1omega+ci*(1-ny**2)*sin3omega))
    Dmat(2,1) =  2.d0*ci*sin1omega*(ci*ny**3*cos1omega + ci*ny*(1-ny**2)*cos3omega + &
                                      nx*nz*(-3.d0*ci*ny**2*sin1omega-ci*(1-ny**2)*sin3omega))
    Dmat(2,3) = -2.d0*ci*sin1omega*(ci*nz**3*cos1omega + ci*nz*(1-nz**2)*cos3omega + &
                                      nx*ny*( 3.d0*ci*nz**2*sin1omega+ci*(1-nz**2)*sin3omega))
    Dmat(3,2) =  2.d0*ci*sin1omega*(ci*nz**3*cos1omega + ci*nz*(1-nz**2)*cos3omega + &
                                      nx*ny*(-3.d0*ci*nz**2*sin1omega-ci*(1-nz**2)*sin3omega))
    Dmat(1,3) =  2.d0*ci*sin1omega*(ci*nx**3*cos1omega + ci*nx*(1-nx**2)*cos3omega + &
                                      ny*nz*(-3.d0*ci*nx**2*sin1omega-ci*(1-nx**2)*sin3omega))
    Dmat(3,1) = -2.d0*ci*sin1omega*(ci*nx**3*cos1omega + ci*nx*(1-nx**2)*cos3omega + &
                                      ny*nz*( 3.d0*ci*nx**2*sin1omega+ci*(1-nx**2)*sin3omega))

    Dmat(1,4) = -ci*sin1omega*(ci*nz*(3.d0-nz**2)*cos1omega + ci*nz*(1+nz**2)*cos3omega - &
                                 4.d0*ci*nx*ny*(ny**2-nx**2)*sin1omega**3)
    Dmat(4,1) =  ci*sin1omega*(ci*nz*(3.d0-nz**2)*cos1omega + ci*nz*(1+nz**2)*cos3omega + &
                                 4.d0*ci*nx*ny*(ny**2-nx**2)*sin1omega**3)
    Dmat(1,5) = -2.d0*sqrt3*sin1omega**2*(nx*ny*(1.d0-nz**2) + nx*ny*(1+nz**2)*cos2omega - &
                                               nz*(ny**2-nx**2)*sin2omega)
    Dmat(5,1) = -2.d0*sqrt3*sin1omega**2*(nx*ny*(1.d0-nz**2) + nx*ny*(1+nz**2)*cos2omega + &
                                               nz*(ny**2-nx**2)*sin2omega)

    Dmat(2,4) = 1.5d0*(nx**2 - ny**2)*ny*nz + ny*nz*(2.d0*ny**2 - 2.d0*nx**2 - 1.d0)*cos2omega + &
                0.5d0*ny*nz*(2.d0*nx**2 + nz**2 + 1.d0)*cos4omega - &
                nx*( 2.d0*ny**2 - nz**2 + (1.d0 - 2.d0*ny**2 + nz**2)*cos2omega)*sin2omega
    Dmat(3,4) = 1.5d0*(nx**2 - ny**2)*nx*nz + nx*nz*(4.d0*ny**2 + 2.d0*nz**2 - 1.d0)*cos2omega - &
                0.5d0*nx*nz*(2.d0*ny**2 + nz**2 + 1.d0)*cos4omega - &
                ny*( 2.d0*nx**2 - nz**2 + (-1.d0 + 2.d0*ny**2 + 3.d0*nz**2)*cos2omega)*sin2omega
    Dmat(4,2) = 1.5d0*(nx**2 - ny**2)*ny*nz + ny*nz*(2.d0*ny**2 - 2.d0*nx**2 - 1.d0)*cos2omega + &
                0.5d0*ny*nz*(2.d0*nx**2 + nz**2 + 1.d0)*cos4omega + &
                nx*( 2.d0*ny**2 - nz**2 + (1.d0 - 2.d0*ny**2 + nz**2)*cos2omega)*sin2omega
    Dmat(4,3) = 1.5d0*(nx**2 - ny**2)*nx*nz + nx*nz*(4.d0*ny**2 + 2.d0*nz**2 - 1.d0)*cos2omega - &
                0.5d0*nx*nz*(2.d0*ny**2 + nz**2 + 1.d0)*cos4omega + &
                ny*( 2.d0*nx**2 - nz**2 + (-1.d0 + 2.d0*ny**2 + 3.d0*nz**2)*cos2omega)*sin2omega

    Dmat(2,5) = 2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(nx*cos1omega - ny*nz*sin1omega)
    Dmat(3,5) =-2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(ny*cos1omega + nx*nz*sin1omega)
    Dmat(5,3) = 2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(ny*cos1omega - nx*nz*sin1omega)
    Dmat(5,2) =-2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(nx*cos1omega + ny*nz*sin1omega)

    Dmat(4,5) = -sqrt3*sin1omega**2*(nx**4 - ny**4 + (nx**2 - ny**2)*(1.d0 + nz**2)*cos2omega - &
                                          4.d0*nx*ny*nz*sin2omega)
    Dmat(5,4) = -sqrt3*sin1omega**2*(nx**4 - ny**4 + (nx**2 - ny**2)*(1.d0 + nz**2)*cos2omega + &
                                          4.d0*nx*ny*nz*sin2omega)
  else if (twoja1 == 7) then 
    Dmat(1,1) = 0.25d0*exp_3omega*(60d0*exp3omega*nx2*ny2*nz2 + 3d0*(1d0-nx2)*(1d0-ny2)*(1d0-nz2) + &
                3d0*exp6omega*(1d0-nx2)*(1d0-ny2)*(1d0-nz2) + &
                5d0*exp2omega*(nx4*(1d0-nx2)+(1d0-nx2)*ny2*nz2 + nx2*(ny4-6d0*ny2*nz2+nz4)) + &
                5d0*exp4omega*(nx4*(1d0-nx2)+(1d0-nx2)*ny2*nz2 + nx2*(ny4-6d0*ny2*nz2+nz4)) + &
                2d0*expomega *(nx6-nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - nx2*(ny4-3d0*ny2*nz2+nz4)) + &
                2d0*exp5omega*(nx6-nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - nx2*(ny4-3d0*ny2*nz2+nz4))   )
    Dmat(2,1) = -2d0*sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*((1d0-nx2)*ny*nz+(1d0+nx2)*ny*nz*cos1omega+nx*(nz2-ny2)*sin1omega)
    Dmat(3,1) = -2d0*sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2*((1d0-ny2)*nx*nz+(1d0+ny2)*nx*nz*cos1omega+ny*(nx2-nz2)*sin1omega)
    Dmat(4,1) = -2d0*sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*((1d0-nz2)*nx*ny+(1d0+nz2)*nx*ny*cos1omega+nz*(ny2-nx2)*sin1omega)
    Dmat(5,1) = 0.25d0*ci*exp_omega_2*(expomega-1d0)*(2d0*nx3*(2d0*nx2+5d0*(1d0-nx2))*cosOmega_2 + &
                                                      nx*(4d0*nx4+5d0*(1d0-nx2)**2)*cos3omega_2 + &
                                                      6d0*nx3*ny2*cos5omega_2 + 3d0*nx*ny4*cos5omega_2 + 6d0*nx3*nz2*cos5omega_2 + &
                                                      6d0*nx*ny2*nz2*cos5omega_2 + 3d0*nx*nz4*cos5omega_2 + &
                                                      30d0*nx2*ny*nz*(nz2-ny2)*sinOmega_2 + &
                                                      5d0*(1d0-3d0*nx2)*ny*nz*(nz2-ny2)*sin3omega_2 + &
                                                      3d0*ny*nz*(ny4-nz4)*sin5omega_2 )
    Dmat(6,1) = -0.5d0*sinOmega_2*(2d0*ny3*(5d0*nx2+2d0*ny2+5d0*nz2)*cosOmega_2 + &
                                            ny*(4d0*ny4+5d0*(1d0-ny2)**2)*cos3omega_2 + &
                                            6d0*ny3*nx2*cos5omega_2 + 3d0*ny*nx4*cos5omega_2 + 6d0*nx2*ny*nz2*cos5omega_2 + &
                                            6d0*ny3*nz2*cos5omega_2 + 3d0*ny*nz4*cos5omega_2 + &
                                            30d0*ny2*nx*nz*(nx2-nz2)*sinOmega_2 + &
                                            5d0*(1d0-3d0*ny2)*nx*nz*(nx2-nz2)*sin3omega_2 + &
                                            3d0*nx*nz*(nz4-nx4)*sin5omega_2 )
    Dmat(7,1) = -0.5d0*sinOmega_2*(2d0*nz3*(2d0*nz2+5d0*(1d0-nz2))*cosOmega_2 + &
                                   nz*(4d0*nz4+5d0*(1d0-nz2)**2)*cos3omega_2 + &
                                   6d0*nz3*nx2*cos5omega_2 + 3d0*nz*nx4*cos5omega_2 + 6d0*nx2*nz*ny2*cos5omega_2 + &
                                   6d0*nz3*ny2*cos5omega_2 + 3d0*nz*ny4*cos5omega_2 + &
                                   30d0*nz2*nx*ny*(nx2-ny2)*sinOmega_2 + &
                                   5d0*(1d0-3d0*nz2)*nx*ny*(ny2-nx2)*sin3omega_2 + &
                                   3d0*nx*ny*(nx4-ny4)*sin5omega_2) 
    Dmat(1,2) = -2d0*sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*((1d0-nx2)*ny*nz+(1d0+nx2)*ny*nz*cos1omega+nx*(ny2-nz2)*sin1omega)
    Dmat(2,2) = 1d0/16d0*exp_3omega*( 3d0*exp2omega*(1d0-5d0*nx2)**2*(1d0-nx2) + 3d0*exp4omega*(1d0-5d0*nx2)**2*(1d0-nx**2) + &
                                       30d0*expomega*nx2*(1d0-nx2)**2 + 30d0*exp5omega*nx2*(1d0-nx2)**2 + &
                                       5d0*(1d0-nx2)**3 + 5d0*exp6omega*(1d0-nx2)**3 + 4d0*exp3omega*nx2*(2d0*nx2-3d0*(1d0-nx2))**2 )
    Dmat(3,2) = -0.25d0*sinOmega_2*(2d0*nz*(6d0*nx4-3d0*nx2*ny2+6d0*ny4+nz4+7d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                    5d0*nz*(9d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cos3omega_2 - 15d0*nx2*ny2*nz*cos5omega_2 + &
                                    5d0*nx2*nz3*cos5omega_2 + 5d0*ny2*nz3*cos5omega_2 + 5d0*nz5*cos5omega_2 + &
                                    2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + &
                                    5d0*nx*ny*(5d0*nx2*ny2-3d0*nz4-3d0*nz2*(1d0-nz2))*sin3omega_2 + &
                                    5d0*nx*ny*(-nx2*ny2+3d0*nz4+3d0*nz2*(1-nz2))*sin5omega_2 )
    Dmat(4,2) = 0.25d0*sinOmega_2*(2d0*ny*(6d0*nx4+nx2*(7d0*ny2-3d0*nz2)+(1d0-nx2)*(ny2+6d0*nz2))*cosOmega_2 + &
                                   5d0*ny*((1-nx2)*ny2+nx2*(ny2+9d0*nz2))*cos3omega_2 - 15d0*nx2*ny*nz2*cos5omega_2 + &
                                   5d0*nx2*ny3*cos5omega_2 + 5d0*ny5*cos5omega_2 + 5d0*ny3*nz2*cos5omega_2 + &
                                   (-2d0)*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(-2d0*nz2+3d0*(1d0-nz2))*sinOmega_2 + &
                                   5d0*nx*nz*(3d0*(1d0-nx2)*ny2+nx2*(3d0*ny2-5d0*nz2))*sin3omega_2 + &
                                   5d0*nx*nz*(-3d0*(1d0-nx2)*ny2 + nx2*(-3d0*ny2+nz2))*sin5omega_2 )
    Dmat(5,2) = -sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*(ny4-nz4+(1d0+nx2)*(ny2-nz2)*cos1omega-4d0*nx*ny*nz*sin1omega)
    Dmat(6,2) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*nz*(2d0*nx4+3d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cosOmega_2 + &
                                                              nz*(-nx2*(ny2-7d0*nz2)+(1d0-nx2)*(2d0*ny2+nz2))*cos3omega_2 + &
                                                              3d0*nx2*ny2*nz*cos5omega_2 + 2d0*ny4*nz*cos5omega_2 + &
                                                            (-1d0)*nx2*nz3*cos5omega_2 + 3d0*ny2*nz3*cos5omega_2 + nz5*cos5omega_2 + &
                                                             2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(nx2-nz2)*sinOmega_2 + &
                                                             nx*ny*(-2d0*ny4-ny2*nz2+nz4+nx2*(3d0*ny2+11d0*nz2))*sin3omega_2 + &
                                                             nx*ny*(2d0*ny4+ny2*nz2-nz4+nx2*(ny2-3d0*nz2))*sin5omega_2 )
    Dmat(7,2) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*ny*(2d0*nx4+nx2*ny2+ny4+(3d0*nx2+ny2)*nz2)*cosOmega_2 + &
                                                              ny*(-nx2*(nz2-7d0*ny2)+(1d0-nx2)*(2d0*nz2+ny2))*cos3omega_2 + &
                                                              3d0*nx2*ny*nz2*cos5omega_2 + 2d0*ny*nz4*cos5omega_2 + &
                                                            (-1d0)*nx2*ny3*cos5omega_2 + 3d0*ny3*nz2*cos5omega_2 + ny5*cos5omega_2 + &
                                                             2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(ny2-nx2)*sinOmega_2 + &
                                                            (-1d0)*nx*nz*((1d0-nx2)*(ny2-2d0*nz2)+nx2*(11d0*ny2+3d0*nz2))*sin3omega_2 + &
                                                             nx*nz*(3d0*nx2*ny2+ny4-2d0*nz4-nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(1,3) = -2d0*sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2*(nx*(1d0-ny2)*nz+nx*(1d0+ny2)*nz*cos1omega+ny*(nz2-nx2)*sin1omega)
    Dmat(2,3) = 0.25d0*sinOmega_2*(2d0*nz*(6d0*nx4-3d0*nx2*ny2+6d0*ny4+nz4+7d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                   5d0*nz*(9d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cos3omega_2 - 15d0*nx2*ny2*nz*cos5omega_2 + &
                                   5d0*nz3*cos5omega_2 - & 
                                   2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + &
                                   5d0*nx*ny*(-5d0*nx2*ny2+3d0*nz4+3d0*nz2*(1d0-nz2))*sin3omega_2 + &
                                   5d0*nx*ny*(nx2*ny2-3d0*nz4-3d0*nz2*(1-nz2))*sin5omega_2 )
    Dmat(3,3) = 1d0/16d0*exp_3omega*(3d0*exp2omega*(1d0-5d0*ny2)**2*(1d0-ny2) + 3d0*exp4omega*(1d0-5d0*ny2)**2*(1d0-ny2) + &
                                     30d0*expomega*ny2*(1d0-ny2)**2 + 30d0*exp5omega*ny2*(1d0-ny2)**2 + 5d0*(1d0-ny2)**3 + &
                                     5d0*exp6omega*(1d0-ny2)**3 + 4d0*exp3omega*ny2*(3d0*nx2-2d0*ny2+3d0*nz2)**2 )
    Dmat(4,3) = 0.125d0*(4d0*ny*nz*(3d0*(3d0*nx4+nx2*(1d0-nx2)-2d0*ny4 + ny2*nz2 - 2d0*nz4) -20d0*ny2*nz2*cos1omega - &
                                    5d0*(3d0*nx4+3d0*nx2*(1d0-nx2)-ny2*nz2)*cos2omega)*sinOmega_2**2 - &
                                    2d0*nx*(nx4+7d0*nx2*(1d0-nx2)+6d0*ny4-33d0*ny2*nz2+6d0*nz4+60d0*ny2*nz2*cos1omega + &
                                    5d0*(nx4+nx2*(1d0-nx2)-3d0*ny2*nz2)*cos2omega) * sin1omega )
    Dmat(5,3) = 1d0/16d0*sqrt15*exp_3omega*(-2d0*(expomega-1d0)**2*(1d0+exp4omega)*nx5*ny + &
                                            2d0*ci*nx4*nz*(-1d0+exp4omega*(2d0*cos2omega-1d0)) + &
                                        8d0*exp3omega*nx3*ny*(-ny2+3d0*nz2+4d0*ny2*cos1omega+(1d0-nx2)*cos2omega)*sinOmega_2**2 - &
                 8d0*exp3omega*nx*ny*(-2d0*ny4+ny2*nz2-3d0*nz4+nz2*(-8d0*ny2*cos1omega + (3d0*ny2+nz2)*cos2omega))*sinOmega_2**2- &
                 4d0*exp3omega*nx2*nz*(5d0*ny2+nz2-4d0*ny2*cos1omega+3d0*(1d0-nx2)*cos2omega)*sin1omega - &
                 4d0*exp3omega*nz*(8d0*ny2*nz2*cos1omega + (ny2-nz2)*(2d0*ny2-nz2-nz2*cos2omega))*sin1omega )
    Dmat(6,3) = sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2 * (nx4-nz4+(1+ny2)*(nx2-nz2)*cos1omega+4d0*nx*ny*nz*sin1omega)
    Dmat(7,3) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)* &
                ( 2d0*nx*cosOmega_2 * ( nx4 - 3d0*nx2*ny2 + 2d0*ny4 + nz2*(1d0+4d0*ny2-nz2) + 4d0*ny2*(2d0*nx2-nz2)*cos1omega + &
                                      (nx2*(nx2-ny2)+2d0*nz4+3d0*nz2*(1d0-nz2))*cos2omega ) + &
                  ny*nz * (2d0*(nx2-ny2)*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + &
                          (nx4 + 3d0*ny2*nz2 - 2d0*nz4 + nx2*(11d0*ny2-nz2))*sin3omega_2 + &
                          (2d0*nz4 + nz2*(1d0-nz2) - nx2*(1d0+2d0*ny2-nz2))*sin5omega_2 ))
    Dmat(1,4) = -2d0*sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*(nx*ny*(1d0-nz2)+nx*ny*(1d0+nz2)*cos1omega+(nx2-ny2)*nz*sin1omega)
    Dmat(2,4) = -0.25d0*sinOmega_2*(2d0*ny*(6d0*nx4+nx2*(7d0*ny2-3d0*nz2) + (1d0-nx2)*(ny2+6d0*nz2))*cosOmega_2 + &
                                    5d0*ny*((1d0-nx2)*ny2+nx2*(ny2+9d0*nz2))*cos3omega_2 + 5d0*nx2*ny3*cos5omega_2 + &
                                    5d0*ny5*cos5omega_2 -15d0*nx2*ny*nz2*cos5omega_2 + 5d0*ny3*nz2*cos5omega_2 + &
                                    2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(-2d0*nz2+3d0*(1d0-nz2))*sinOmega_2 - &
                                    5d0*nx*nz*(3d0*(1d0-nx2)*ny2+nx2*(3d0*ny2-5d0*nz2))*sin3omega_2 + &
                                    5d0*nx*nz*(3d0*(1d0-nx2)*ny2+nx2*(3d0*ny2-nz2))*sin5omega_2 )
    Dmat(3,4) = 0.125d0*(4d0*ny*nz*(3d0*(3d0*nx4+nx2*(1d0-nx2)-2d0*ny4+ny2*nz2-2d0*nz4) - 20d0*ny2*nz2*cos1omega - &
                                     5d0*(3d0*nx4+3d0*nx2*(1d0-nx2)-ny2*nz2)*cos2omega)*sinOmega_2**2 + &
                          2d0*nx*(nx4+7d0*nx2*(1d0-nx2)+6d0*ny4-33d0*ny2*nz2+6d0*nz4+60*ny2*nz2*cos1omega + &
                          5d0*(nx4+nx2*(1d0-nx2)-3d0*ny2*nz2)*cos2omega)*sin1omega )
    Dmat(4,4) = 1d0/16d0*exp_3omega*(3d0*exp2omega*(1d0-5d0*nz2)**2*(1d0-nz2)+3d0*exp4omega*(1d0-5d0*nz2)**2*(1d0-nz2) + &
                                     30d0*expomega*nz2*(1d0-nz2)**2+30d0*exp5omega*nz2*(1d0-nz2)**2 + &
                                     5d0*(1d0+exp6omega)*(1d0-nz2)**3 + 4d0*exp3omega*(-2d0*nz3+3d0*nz*(1d0-nz2))**2) 
    Dmat(5,4) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*ny*cosOmega_2*(ny4-3d0*ny2*nz2+2d0*nz4+nx2*(ny2+5d0*nz2) - &
                                        4d0*nz2*(1d0-3d0*ny2-nz2)*cos1omega + (2d0*nx4+3d0*nx2*(1d0-nx2)+ny2*(ny2-nz2))*cos2omega)+&
                                        nx*nz*(2d0*(ny2-nz2)*(3d0-5d0*nz2)*sinOmega_2 + &
                                        (-2d0*nx4+ny4+11d0*ny2*nz2-nx2*(ny2-3d0*nz2))*sin3omega_2 + &
                                        (2d0*nx4+nx2*(1d0-nx2)-ny2*(ny2+3d0*nz2))*sin5omega_2))
    Dmat(6,4) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*nx*cosOmega_2*(nx4-3d0*nx2*nz2+2d0*nz4+ny2*(nx2+5d0*nz2) + &
                                        4d0*nz2*(2d0*nx2-ny2)*cos1omega + (nx4+2d0*ny4+3d0*ny2*nz2+nx2*(3d0*ny2-nz2))*cos2omega)+&
                                        ny*nz*(-2d0*(nx2-nz2)*(3d0-5d0*nz2)*sinOmega_2 - &
                                        (nx4-nx2*ny2-2d0*ny4+(11d0*nx2+3d0*ny2)*nz2)*sin3omega_2 + &
                                        (nx4-nx2*ny2-2d0*ny4+nz2*(3d0*nx2-ny2))*sin5omega_2))
    Dmat(7,4) = -sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*(nx4-ny4+(nx2-ny2)*(1d0+nz2)*cos1omega-4d0*nx*ny*nz*sin1omega)
    Dmat(1,5) = -0.25d0*ci*exp_omega_2*(expomega-1d0)*(2d0*nx3*(2d0*nx2+5d0*(1d0-nx2))*cosOmega_2 + &
                                                      nx*(4d0*nx4+5d0*(1d0-nx2)**2)*cos3omega_2 + &
                                                      6d0*nx3*ny2*cos5omega_2 + 3d0*nx*ny4*cos5omega_2 + 6d0*nx3*nz2*cos5omega_2 + &
                                                      6d0*nx*ny2*nz2*cos5omega_2 + 3d0*nx*nz4*cos5omega_2 + &
                                                      40d0*nx2*ny*nz*(ny2-nz2)*sinOmega_2**3 + &
                                                      8d0*nz*(ny5-ny*nz4)*(2d0+3d0*cos1omega)*sinOmega_2**3)
    Dmat(2,5) = -sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*(ny4-nz4+(1d0+nx2)*(ny2-nz2)*cos1omega+4d0*nx*ny*nz*sin1omega)
    Dmat(3,5) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*&
        (-2d0*nz*cosOmega_2*(2d0*ny4-3d0*ny2*nz2+nz4+nx2*(5d0*ny2+nz2)-4d0*ny2*(1d0-ny2-3d0*nz2)*cos1omega + &
         (2d0*nx4+3d0*nx2*(1d0-nx2)-ny2*nz2+nz4)*cos2omega) + &
         nx*ny*(2d0*(ny2-nz2)*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + (2d0*nx4+nx2*(-3d0*ny2+nz2)-nz2*(11d0*ny2+nz2))*sin3omega_2 + &
         (-2d0*nx4-nx2*(1d0-nx2)+3d0*ny2*nz2+nz4)*sin5omega_2 ))
    Dmat(4,5) = 1d0/16d0*sqrt15*exp_3omega*( 2d0*(expomega-1d0)**2*(1d0+exp4omega)*nx5*nz - &
                                              2d0*ci*nx4*ny*(exp4omega*(2d0*cos2omega-1d0)-1d0) - &
                                    8d0*exp3omega*nx3*nz*(3d0*ny2-nz2+4d0*nz2*cos1omega + (1d0-nx2)*cos2omega)*sinOmega_2**2 + &
                8d0*exp3omega*nx*nz*(-3d0*ny4+ny2*nz2-2d0*nz4+ny2*(-8d0*nz2*cos1omega+(ny2+3d0*nz2)*cos2omega))*sinOmega_2**2 + &
                4d0*exp3omega*nx2*ny*(ny2+5d0*nz2-4d0*nz2*cos1omega+3d0*(1d0-nx2)*cos2omega)*sin1omega + &
                4d0*exp3omega*ny*(8d0*ny2*nz2*cos1omega+(ny2-nz2)*(ny2-2d0*nz2+ny2*cos2omega))*sin1omega )
    Dmat(5,5) = 1d0/16d0*exp_3omega*(60d0*exp3omega*nx2*(ny2-nz2)**2+3d0*(1d0-nx2)*(4d0*nx4+4d0*nx2*(1d0-nx2)+(ny2-nz2)**2) + &
                                    3d0*exp6omega*(1d0-nx2)*(4d0*nx4+4d0*nx2*(1d0-nx2)+(ny2-nz2)**2) + &
                            2d0*expomega*((-2d0*nx3+nx*ny2)**2-2d0*(2d0*nx4+9d0*nx2*ny2-8d0*ny4)*nz2 + nz4*(1d0+15d0*ny2-nz2)) + &
                            2d0*exp5omega*((-2d0*nx3+nx*ny2)**2-2d0*(2d0*nx4+9d0*nx2*ny2-8d0*ny4)*nz2+ nz4*(1d0+15d0*ny2-nz2)) + &
                            5d0*exp2omega*(4d0*nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - 4d0*nx2*(ny4-6d0*ny2*nz2+nz4)) + &
                            5d0*exp4omega*(4d0*nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - 4d0*nx2*(ny4-6d0*ny2*nz2+nz4)))
    Dmat(6,5) = -0.25d0*sinOmega_2*( 2*nz*(-15d0*nx2*ny2+nz4-5d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                    nz*(-5d0*(2d0*nx4-5d0*nx2*ny2+2d0*ny4)-3d0*nz4+5d0*nz2*(1d0-nz2))*cos3omega_2 + &
                                    6d0*nx4*nz*cos5omega_2 -3d0*nx2*ny2*nz*cos5omega_2 + 6d0*ny4*nz*cos5omega_2 - &
                                    3d0*nx2*nz3*cos5omega_2 - &
                                    3d0*ny2*nz3*cos5omega_2 -3d0*nz5*cos5omega_2 + 30d0*nx*ny*(nx2-nz2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*ny*(2d0*nx4+nx2*ny2+2*ny4-3d0*nz4+5d0*nz2*(1d0-nz2))*sin3omega_2 + &
                                    3d0*nx*ny*(2d0*nx4+5d0*nx2*ny2+2d0*ny4+5d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(7,5) = -0.25d0*sinOmega_2*( 2d0*ny*(-ny4+5d0*ny2*nz2+5d0*nx2*(ny2+3d0*nz2))*cosOmega_2 + &
                                    ny*(10d0*nx4-5d0*nx2*ny2+3d0*ny4-5d0*(5d0*nx2+ny2)*nz2+10d0*nz4)*cos3omega_2 - &
                                    6d0*nx4*ny*cos5omega_2 + 3d0*nx2*ny*nz2*cos5omega_2 - 6d0*ny*nz4*cos5omega_2 + &
                                    3d0*nx2*ny3*cos5omega_2 + &
                                    3d0*ny3*nz2*cos5omega_2 + 3d0*ny5*cos5omega_2 - 30d0*nx*nz*(nx2-ny2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2*nz4-3d0*ny4+nz2*(1d0+4d0*ny2-nz2))*sin3omega_2 + &
                                    3d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2d0*nz4+5d0*ny4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(1,6) = 0.5d0*sinOmega_2*(2d0*ny3*(5d0*nx2+2d0*ny2+5d0*nz2)*cosOmega_2 + &
                                            ny*(4d0*ny4+5d0*(1d0-ny2)**2)*cos3omega_2 + &
                                            6d0*ny3*nx2*cos5omega_2 + 3d0*ny*nx4*cos5omega_2 + 6d0*nx2*ny*nz2*cos5omega_2 + &
                                            6d0*ny3*nz2*cos5omega_2 + 3d0*ny*nz4*cos5omega_2 + &
                                            30d0*ny2*nx*nz*(nz2-nx2)*sinOmega_2 + &
                                            5d0*(1d0-3d0*ny2)*nx*nz*(nz2-nx2)*sin3omega_2 + &
                                            3d0*nx*nz*(nx4-nz4)*sin5omega_2 )
    Dmat(2,6) = 0.25d0*sqrt15*sinOmega_2*( 2d0*nz*(2d0*nx4+3d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cosOmega_2 + &
                                           nz*(-nx2*(ny2-7d0*nz2)+(1d0-nx2)*(2d0*ny2+nz2))*cos3omega_2 + &
                                           3d0*nx2*ny2*nz*cos5omega_2 + 2d0*ny4*nz*cos5omega_2 - nx2*nz3*cos5omega_2 + &
                                           3d0*ny2*nz3*cos5omega_2 + nz5*cos5omega_2 - &
                                           2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(nx2-nz2)*sinOmega_2 - &
                                           nx*ny*(-2d0*ny4-ny2*nz2+nz4+nx2*(3d0*ny2+11d0*nz2))*sin3omega_2 + &
                                           nx*ny*(-2d0*ny4-ny2*nz2+nz4-nx2*(ny2-3d0*nz2))*sin5omega_2 )
    Dmat(3,6) = sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2*(nx4-nz4+(1d0+ny2)*(nx2-nz2)*cos1omega-4d0*nx*ny*nz*sin1omega)
    Dmat(4,6) =-0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*nx*cosOmega_2*(nx4-3d0*nx2*nz2+2d0*nz4+ny2*(nx2+5d0*nz2) + &
                                        4d0*nz2*(2d0*nx2-ny2)*cos1omega + (nx4+2d0*ny4+3d0*ny2*nz2+nx2*(3d0*ny2-nz2))*cos2omega)+&
                                        ny*nz*(2d0*(nx2-nz2)*(3d0-5d0*nz2)*sinOmega_2 + &
                                        (nx4-nx2*ny2-2d0*ny4+(11d0*nx2+3d0*ny2)*nz2)*sin3omega_2 - &
                                        (nx4-nx2*ny2-2d0*ny4+nz2*(3d0*nx2-ny2))*sin5omega_2))
    Dmat(5,6) = -0.25d0*sinOmega_2*(-2*nz*(-15d0*nx2*ny2+nz4-5d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                    nz*(5d0*(2d0*nx4-5d0*nx2*ny2+2d0*ny4)+3d0*nz4-5d0*nz2*(1d0-nz2))*cos3omega_2 - &
                                    6d0*nx4*nz*cos5omega_2 +3d0*nx2*ny2*nz*cos5omega_2 - 6d0*ny4*nz*cos5omega_2 + &
                                    3d0*nx2*nz3*cos5omega_2 + &
                                    3d0*ny2*nz3*cos5omega_2 +3d0*nz5*cos5omega_2 + 30d0*nx*ny*(nx2-nz2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*ny*(2d0*nx4+nx2*ny2+2*ny4-3d0*nz4+5d0*nz4*(1d0-nz2))*sin3omega_2 + &
                                    3d0*nx*ny*(2d0*nx4+5d0*nx2*ny2+2d0*ny4+5d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(6,6) = 1d0/16d0*exp_3omega*(60d0*exp3omega*ny2*(1d0-ny2-2d0*nz2)**2 + &
                                     3d0*(1d0-ny2)*(nx4+nx2*(4d0*ny2-2d0*nz2)+(2d0*ny2+nz2)**2) + &
                                     3d0*exp6omega*(1d0-ny2)*(nx4+nx2*(4d0*ny2-2d0*nz2)+(2d0*ny2+nz2)**2) + &
                                2d0*expomega*(nx4*(ny2+16d0*nz2)+(-2d0*ny3+ny*nz2)**2-2d0*nx2*(2d0*ny4+9d0*ny2*nz2-8d0*nz4)) + &
                                2d0*exp5omega*(nx4*(ny2+16d0*nz2)+(-2d0*ny3+ny*nz2)**2-2d0*nx2*(2d0*ny4+9d0*ny2*nz2-8d0*nz4)) + &
                    5d0*(exp2omega+exp4omega)*(nx6-nx4*(4d0*ny2+nz2)+(-2d0*ny2*nz+nz3)**2+nx2*(4d0*ny4+24d0*ny2*nz2-nz4)) )
    Dmat(7,6) = -0.25d0*sinOmega_2*(2d0*nx*(nx4-5d0*nx2*(1d0-nx2)-15d0*ny2*nz2)*cosOmega_2 + &
                                    nx*(-3d0*nx4+5d0*nx2*(1d0-nx2)-5d0*(2d0*ny4-5d0*ny2*nz2+2d0*nz4))*cos3omega_2 + &
                                    cos5omega_2*(-3d0*nx5-3d0*nx3*ny2+6d0*nx*ny4-3d0*nx3*nz2-3d0*nx*ny2*nz2+6d0*nx*nz4) + &
                                    30d0*ny*nz*(nx2-ny2)*(nx2-nz2)*sinOmega_2 + &
                                    5d0*ny*nz*(-3d0*nx4+5d0*nx2*(1d0-nx2)+2d0*ny4+ny2*nz2+2d0*nz4)*sin3omega_2 + &
                                    3d0*ny*nz*(5d0*nx4+5d0*nx2*ny2+2d0*ny4+2d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(1,7) =  0.5d0*sinOmega_2*(2d0*nz3*(2d0*nz2+5d0*(1d0-nz2))*cosOmega_2 + &
                                   nz*(4d0*nz4+5d0*(1d0-nz2)**2)*cos3omega_2 + &
                                   6d0*nz3*nx2*cos5omega_2 + 3d0*nz*nx4*cos5omega_2 + 6d0*nx2*nz*ny2*cos5omega_2 + &
                                   6d0*nz3*ny2*cos5omega_2 + 3d0*nz*ny4*cos5omega_2 + &
                                   30d0*nz2*nx*ny*(nx2-ny2)*sinOmega_2 + &
                                   5d0*(1d0-3d0*nz2)*nx*ny*(nx2-ny2)*sin3omega_2 + &
                                   3d0*nx*ny*(ny4-nx4)*sin5omega_2) 
    Dmat(2,7) =-0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*ny*(2d0*nx4+nx2*ny2+ny4+(3d0*nx2+ny2)*nz2)*cosOmega_2 + &
                                                              ny*(-nx2*(nz2-7d0*ny2)+(1d0-nx2)*(2d0*nz2+ny2))*cos3omega_2 + &
                                                              3d0*nx2*ny*nz2*cos5omega_2 + 2d0*ny*nz4*cos5omega_2 + &
                                                            (-1d0)*nx2*ny3*cos5omega_2 + 3d0*ny3*nz2*cos5omega_2 + ny5*cos5omega_2 + &
                                                             2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(nx2-ny2)*sinOmega_2 + &
                                                            nx*nz*((1d0-nx2)*(ny2-2d0*nz2)+nx2*(11d0*ny2+3d0*nz2))*sin3omega_2 + &
                                                             nx*nz*(-3d0*nx2*ny2-ny4+2d0*nz4+nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(3,7) = 1d0/16d0*sqrt15*exp_3omega*(-ci*(exp2omega-1d0)*(exp2omega+1d0)**2*nx5 - &
                                            (expomega-1d0)**2*(1d0-6d0*exp2omega+exp4omega)*nx4*ny*nz + &
                                8d0*exp3omega*nx2*ny*nz*(ny2-3d0*nz2-8d0*ny2*cos1omega+(3d0*ny2-nz2)*cos2omega)*sinOmega_2**2 - &
                            8d0*exp3omega*ny*nz*(ny2*(2d0*ny2-nz2+4d0*nz2*cos1omega)+nz2*(ny2+2d0*nz2)*cos2omega)*sinOmega_2**2 + &
                            4d0*exp3omega*nx3*(-3d0*ny2+nz2+8d0*ny2*cos1omega-(ny2-3d0*nz2)*cos2omega)*sin1omega + &
                            4d0*exp3omega*nx*(ny2*(2d0*ny2+5d0*nz2-4d0*nz2*cos1omega)+nz2*(3d0*ny2+2d0*nz2)*cos2omega)*sin1omega )
    Dmat(4,7) = -sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*(nx4-ny4+(nx2-ny2)*(1d0+nz2)*cos1omega+4d0*nx*ny*nz*sin1omega)
    Dmat(5,7) = -0.25d0*sinOmega_2*(-2d0*ny*(-ny4+5d0*ny2*nz2+5d0*nx2*(ny2+3d0*nz2))*cosOmega_2 - &
                                    ny*(10d0*nx4-5d0*nx2*ny2+3d0*ny4-5d0*(5d0*nx2+ny2)*nz2+10d0*nz4)*cos3omega_2 + &
                                    6d0*nx4*ny*cos5omega_2 - 3d0*nx2*ny*nz2*cos5omega_2 + 6d0*ny*nz4*cos5omega_2 - &
                                    3d0*nx2*ny3*cos5omega_2 - &
                                    3d0*ny3*nz2*cos5omega_2 - 3d0*ny5*cos5omega_2 - 30d0*nx*nz*(nx2-ny2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2*nz4-3d0*ny4+nz2*(1d0+4d0*ny2-nz2))*sin3omega_2 + &
                                    3d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2d0*nz4+5d0*ny4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(6,7) = -0.25d0*sinOmega_2*(-2d0*nx*(nx4-5d0*nx2*(1d0-nx2)-15d0*ny2*nz2)*cosOmega_2 + &
                                    nx*(3d0*nx4-5d0*nx2*(1d0-nx2)+5d0*(2d0*ny4-5d0*ny2*nz2+2d0*nz4))*cos3omega_2 + &
                                    cos5omega_2*(3d0*nx5+3d0*nx3*ny2-6d0*nx*ny4+3d0*nx3*nz2+3d0*nx*ny2*nz2-6d0*nx*nz4) + &
                                    30d0*ny*nz*(nx2-ny2)*(nx2-nz2)*sinOmega_2 + &
                                    5d0*ny*nz*(-3d0*nx4+5d0*nx2*(1d0-nx2)+2d0*ny4+ny2*nz2+2d0*nz4)*sin3omega_2 + &
                                    3d0*ny*nz*(5d0*nx4+5d0*nx2*ny2+2d0*ny4+2d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(7,7) = 1d0/16d0*exp_3omega*(60d0*exp3omega*nz2*(1d0-2d0*ny2-nz2)**2 + &
                                     2d0*expomega*((nx4-18d0*nx2*ny2+ny4)*nz2+4d0*nz6+16d0*nx2*ny2*(1d0-nz2)-4d0*nz4*(1d0-nz2)) + &
                                     2d0*exp5omega*((nx4-18d0*nx2*ny2+ny4)*nz2+4d0*nz6+16d0*nx2*ny2*(1d0-nz2)-4d0*nz4*(1d0-nz2)) + &
                                     3d0*(1d0+exp6omega)*(1d0-nz2)*(4d0*nz4+4d0*nz2*(1d0-nz2)+(1d0-2d0*ny2-nz2)**2) + &
                        5d0*(exp2omega+exp4omega)*(-4d0*(nx4-6d0*nx2*ny2+ny4)*nz2+4d0*nz4*(1d0-nz2)+(1d0-nz2)*(1d0-2d0*ny2-nz2)**2) )

  else
    write(*,*) " This subroutine can only deal with J = 1/2, 1, 2, 3 !!! "
  end if

    ! using the same order of orbitals as Wannier90 
    ! l=0 : s
    ! l=1 : pz, px, py
    ! l=2 : dz2, dxz, dyz, dx2-y2, dxy
  !if (twoja1 == 3) then
  !  Dtmp(1,:) = Dmat(1,:)
  !  Dmat(1,:) = Dmat(3,:)
  !  Dmat(3,:) = Dtmp(1,:)
  !  Dtmp(:,1) = Dmat(:,1)
  !  Dmat(:,1) = Dmat(:,3)
  !  Dmat(:,3) = Dtmp(:,1)
  !elseif (twoja1 == 5) then 
  !  Dtmp(1,:) = Dmat(1,:)
  !  Dmat(1,:) = Dmat(5,:)
  !  Dmat(5,:) = Dtmp(1,:)
  !  Dtmp(:,1) = Dmat(:,1)
  !  Dmat(:,1) = Dmat(:,5)
  !  Dmat(:,5) = Dtmp(:,1)
  !  Dtmp(2,:) = Dmat(2,:)
  !  Dmat(2,:) = Dmat(3,:)
  !  Dmat(3,:) = Dtmp(2,:)
  !  Dtmp(:,2) = Dmat(:,2)
  !  Dmat(:,2) = Dmat(:,3)
  !  Dmat(:,3) = Dtmp(:,2)
  !endif 

end subroutine Dmatrix

subroutine O3_to_axis_angle(R, axis, angle, detR)
  use iso_fortran_env, only: dp=>real64
  implicit none
  real(dp), intent(in)  :: R(3,3)
  real(dp), intent(out) :: axis(3), angle, detR

  real(dp), parameter :: ONE=1.0_dp, ZERO=0.0_dp, TWO=2.0_dp
  real(dp), parameter :: eps=1.0e-4_dp, eps_ang=1.0e-4_dp
  real(dp), parameter :: pi = acos(-ONE)
  real(dp) :: Rrot(3,3), tr, cosang, sinang, r11,r22,r33, nrm

  detR = R(1,1)*(R(2,2)*R(3,3) - R(2,3)*R(3,2)) &
       - R(1,2)*(R(2,1)*R(3,3) - R(2,3)*R(3,1)) &
       + R(1,3)*(R(2,1)*R(3,2) - R(2,2)*R(3,1))

  if (detR < 0.0_dp) then
    Rrot = -R
  else
    Rrot =  R
  end if

  tr = Rrot(1,1) + Rrot(2,2) + Rrot(3,3)
  cosang = (tr - ONE)/TWO
  if (cosang >  ONE) cosang =  ONE
  if (cosang < -ONE) cosang = -ONE
  angle = acos(cosang)

  if (abs(angle) < eps_ang) then
    axis = (/ZERO, ZERO, ONE/)
    return
  end if

  if (abs(angle - pi) < eps_ang) then
    r11=Rrot(1,1); r22=Rrot(2,2); r33=Rrot(3,3)
    if (r11 >= r22 .and. r11 >= r33) then
      axis(1)=sqrt(max(ZERO,(r11+ONE)/TWO))
      axis(2)=(Rrot(1,2)+Rrot(2,1))/max(eps,4.0_dp*axis(1))
      axis(3)=(Rrot(1,3)+Rrot(3,1))/max(eps,4.0_dp*axis(1))
    else if (r22 >= r11 .and. r22 >= r33) then
      axis(2)=sqrt(max(ZERO,(r22+ONE)/TWO))
      axis(1)=(Rrot(1,2)+Rrot(2,1))/max(eps,4.0_dp*axis(2))
      axis(3)=(Rrot(2,3)+Rrot(3,2))/max(eps,4.0_dp*axis(2))
    else
      axis(3)=sqrt(max(ZERO,(r33+ONE)/TWO))
      axis(1)=(Rrot(1,3)+Rrot(3,1))/max(eps,4.0_dp*axis(3))
      axis(2)=(Rrot(2,3)+Rrot(3,2))/max(eps,4.0_dp*axis(3))
    end if
  else

    sinang = sin(angle)
    axis(1) = (Rrot(3,2) - Rrot(2,3)) / (TWO*sinang)
    axis(2) = (Rrot(1,3) - Rrot(3,1)) / (TWO*sinang)
    axis(3) = (Rrot(2,1) - Rrot(1,2)) / (TWO*sinang)
  end if

  nrm = sqrt(sum(axis*axis))
  if (nrm > eps) axis = axis/nrm

end subroutine O3_to_axis_angle



  subroutine cross(a, b, c)
    implicit none
    real(dp), intent(in) :: a(3), b(3)
    real(dp), intent(out) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end subroutine cross

endmodule