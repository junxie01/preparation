! 2017/06/04
! id1 is the id of signal
! id2 is the id of frequency
subroutine filter(sigm,id1,id2)
use globe_data,only : plan,nn,FFTW_FORWARD,FFTW_ESTIMATE,npts,nk,czero,plan2,FFTW_BACKWARD,dt,nmax,sig
implicit none
integer::    id1,id2,itest
real(4)::    sigm(nmax)
complex(8):: s(nmax),sf(nmax)
s=czero
sf=czero
!            do itest=1,npts
!               write(*,*)sig(itest,1)
!            enddo
!write(*,*)id1,id2,npts
s(1:npts) = dcmplx(dble(sig(1:npts,id1)),0d0)
!write(*,*)' do fft'
!write(*,*)plan,nn,FFTW_FORWARD,FFTW_ESTIMATE
call dfftw_plan_dft_1d(plan,nn,s,sf,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
!do itest=1,nk
!   write(*,*)sf(itest)
!enddo
!write(*,*)nn,nk,dom

sf(1) = sf(1)/2.0d0
sf(nk) = dcmplx(dreal(sf(npts)),0.0d0)
sf(1:nk)=sf(1:nk)*dble(dt)
call taperf(sf,id2)
s=czero
sf(nk+1:nn)=czero
!write(*,*)' do ifft'
! make forward FFT for seismogram: sf ==> s
call dfftw_plan_dft_1d(plan2,nn,sf,s,FFTW_BACKWARD, FFTW_ESTIMATE)
call dfftw_execute(plan2)
call dfftw_destroy_plan(plan2)
sigm(1:npts)=2.0*real(dreal(s(1:npts)))/nn/dt
return
end subroutine
