   subroutine pickconfig(i,r,N,T,pos,Vtot)
   implicit none
   integer, intent(in) :: i,N
   real*8, intent(in) :: T,pos(6,3),Vtot
   real*8, intent(inout) :: r(N,3)
   real*8 Z,Vtemp(6),dtemp,rand,prob
   integer j,k,count,pickk
   
   ! calculate potential for each possible configuration
   ! for the ith atom\
   do k=1,6
      Vtemp(k) = 0
      do j=1,i-1
         dtemp = 0
         dtemp = (r(j,1)-pos(k,1))**2 + (r(j,2)-pos(k,2))**2 + (r(j,3)-pos(k,3))**2
	 print *, dtemp,j,k
         Vtemp(k) = Vtemp(k) + 4*((1.d0/dtemp)**6 - (1.d0/dtemp)**3)
      enddo 
      Vtemp(k) = Vtemp(k) + Vtot
      print *, "Vtemp = ", Vtemp(k),k
   enddo 
   
   ! calculate partition function
   Z = 0
   do k=1,6
      Z = Z + exp(-Vtemp(k)/T)
   enddo 
   print *, "Z = ", Z
   ! check to make sure Z isn't zero
   ! if zero, try new configuration
   !do while (Z==0)
   !   call pickdir(i,r,N,pos)
      
   !enddo  
   
   ! generate random number 
   call random_number(rand)
   ! pick new position based on random number
   ! Metropolis algorithm
   count = 1
   prob = 0
   do while (count < 7)
      prob = prob + exp(-Vtemp(count)/T)/Z
      !print *, "prob = ", prob, rand, exp(-Vtemp(count)/T)/Z, count
      if (rand <= prob) then
         !print *, rand, exp(-Vtemp(count)/T)/Z, count
         pickk = count
	 count = 20
      else
         count = count + 1
      end if 
   enddo 
   
   ! define the new position
   r(i,1) = pos(pickk,1)
   r(i,2) = pos(pickk,2)
   r(i,3) = pos(pickk,3)   
   
   end subroutine pickconfig