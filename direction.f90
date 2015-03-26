! choose a direction for the next atom
! randomly choose theta and phi

subroutine pickdir(i,r,N,pos)
   implicit none
   integer, intent(in) :: i,N
   real*8, intent(inout) :: r(N,3),pos(6,3)
   real*8 theta,phi,pi
   integer j
   
   ! define pi
   pi = 4.d0 * datan(1.d0)
   
   !random number
   call random_number(theta)
   call random_number(phi)
   !print *, theta,phi
   
   ! sample theta from cos(theta)
   ! theta between 0 and pi
   ! phi between 0 and 2pi
   theta = 2*theta - 1 ! between 1 and -1
   theta = acos(theta) ! becomes an angle
   phi = phi*2.d0*pi
   !print *, theta, phi
   
   ! calculate possible new positions, with respect to zero
   ! define theta
   ! define 6 possible phi
   do j=1,6
     phi = (j-1)*pi/6.d0 + phi
     pos(j,1) = r(i-1,1)+sin(theta)*cos(phi)
     pos(j,2) = r(i-1,2)+sin(theta)*sin(phi)
     pos(j,3) = r(i-1,3)+cos(theta)
     !print *, pos(j,1),pos(j,2),pos(j,3)
   enddo 
   
   ! calcluate the energy of each of the new positions
   !call energycalc(i,r,N,pos,T)
   
end subroutine pickdir
