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
