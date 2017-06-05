! 2017/06/04
subroutine do_norm3
use globe_data
implicit none
integer i
real temp
real x1(nmax),x2(nmax),x3(nmax)
real tx1(nmax),tx2(nmax),tx3(nmax)
!halfwidth=int(1.0/fre1(2)/dt/2.0) ! Besen et al. (2007) half of the maximum period of the passband filter
!halfwidth=20                    ! from Yingjie Yang's code
!halfwidth=int(64/dt)
call filter(x1,1,1) ! filter the waveform
call filter(x2,2,1)
call filter(x3,3,1)
!write(*,*)'Smooth the waveform'
call smooth(x1) ! get the running absolute average waveform
call smooth(x2)
call smooth(x3)
call filter(tx1,1,2) ! filter the waveform
call filter(tx2,2,2)
call filter(tx3,3,2)
do i=1,npts 
   temp=max(x1(i),x2(i),x3(i))
   if(comb.eq.1)then
      sig(i,1)=tx1(i)/temp
      sig(i,2)=tx2(i)/temp
      sig(i,3)=tx3(i)/temp
   else
      sig(i,1)=tx1(i)/x1(i)
      sig(i,2)=tx2(i)/x2(i)
      sig(i,3)=tx3(i)/x3(i)
   endif
enddo
end subroutine
