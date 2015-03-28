! main program for advanced monte carlo
! creates a chain polymer
! in 3 dimensions

program chainmc
   implicit none
   integer N,i
   real*8 T,pos,Vtot,maxx,maxy,maxz,minx,miny,minz,diffx,diffy,diffz
   real*8, allocatable, dimension(:,:) :: r
   real*8, allocatable, dimension(:,:) :: V
   real*8, allocatable, dimension(:,:) :: temppos
   
   ! number atoms in the chain
   N=100
   
   ! Temperature
   T = 2.d0
   
   ! seed random number generator
   call random_seed()
   
   ! allocate memory for the positions in the chain
   allocate (r(N,3))
   ! allocate memory for the potentials between pairs
   allocate (V(N,N))
   ! initilize to zero
   r=0.d0
   V=0.d0
   
   ! set positions of the first two atoms (0,0,0) and (1,0,0)
   ! perhaps put this in a subroutine to keep clean
   r(1,1)=0
   r(1,2)=0
   r(1,3)=0
   r(2,1)=1
   r(2,2)=0
   r(2,3)=0
   
   ! add the rest of the atoms
   do i=3,N
      print *,i
      ! calculate the energy of the previous configuration
      call calcpot(i,V,Vtot,N,r)
      
      ! pick new direction, this will fill positions (r)
      ! first pick possible directions
      allocate (temppos(6,3))
      call pickdir(i,r,N,temppos)
      
      ! then calculate the energy of each configuration
      ! use metropolis algorithm to decide "best" configuration
      call pickconfig(i,r,N,T,temppos,Vtot)
      deallocate(temppos)
   enddo
   
   ! check the lengths of all adjacent segements
   do i=1,N-1
      pos = 0
      pos = (r(i,1)-r(i+1,1))**2 + (r(i,2)-r(i+1,2))**2 + (r(i,3)-r(i+1,3))**2
      pos = sqrt(pos)
      !if (pos > 1) then
      !   print *, i,i+1, " not unit length"
      !end if 
      ! print *, i,i+1,pos
   enddo 
   
   !print all positions
   do i=1,N
      print *, r(i,1),r(i,2),r(i,3)
   enddo
   
   ! calculate the size of the system
   ! should be a subroutine
   call calcsize(r,N)
   
   ! prepare to run program many times to get many chains
   ! for some sort of error
   deallocate(r)
   deallocate(V)


end program chainmc
