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
