! 2017/06/04
! one component normilization
subroutine do_norm1
use globe_data
use sacio
implicit none
type(sac_head) :: sacheada
integer i,nerr,itest
real(4),dimension(nmax):: x1,tx1
!halfwidth=int(1.0/fre1(2)/dt/2.0) !Besen et al. (2007) half of the maximum period of the passband filter
!halfwidth=20                  ! from Yingjie Yang's code
!halfwidth=int(64/dt)          ! Fanchi Lin 2007 
!            do itest=1,npts
!               write(*,*)sig(itest,1)
!            enddo
call filter(x1,1,1)            ! filter the waveform
!do i=1,npts
!   write(*,*)x1(i)
!enddo
call initial_sachead(sacheada)
sacheada%delta=dt
sacheada%npts=npts
!call write_sac('test.bp',x1,sacheada,nerr)
!stop
call smooth(x1)      ! get the running absolute average waveform
call filter(tx1,1,2)           ! filter the waveform
do i=1,npts 
   sig(i,1)=tx1(i)/x1(i)
enddo
end subroutine
