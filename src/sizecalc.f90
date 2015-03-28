   subroutine calcsize(r,N)
   implicit none
   integer, intent(in) :: N
   real*8, intent(in) :: r(N,3)
   real*8 maxx,maxy,maxz,minx,miny,minz,diffx,diffy,diffz
   real*8 rsqvar,rvar,rstd,rmm
   integer j
   
   ! max/min value for each coordinate
   maxx = maxval(r(:,1))
   minx = minval(r(:,1))
   maxy = maxval(r(:,2))
   miny = minval(r(:,2))
   maxz = maxval(r(:,3))
   minz = minval(r(:,3))
   
   ! calculate the distance in each dimension
   diffx = maxx - minx
   diffy = maxy - miny
   diffz = maxz - minz
   rmm = 0.5*sqrt(diffx**2 + diffy**2 + diffz**2)
   
   ! std. deviation?
   rsqvar = 0
   rvar = 0
   do j=1,N
      rsqvar = rsqvar + (r(j,1)**2 + r(j,2)**2 + r(j,3)**2)
      rvar = rvar + sqrt(r(j,1)**2 + r(j,2)**2 + r(j,3)**2)
   enddo
   rsqvar = rsqvar/real(N)
   rvar = rvar/real(N)
   rstd = sqrt(rsqvar - rvar**2)
   print *, rmm, rsqvar, rvar, rstd
   print *, log(rmm),log(rstd),log(real(N))
      
   end subroutine calcsize