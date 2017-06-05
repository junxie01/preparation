!2017/06/04
subroutine do_whiten2
use globe_data,only: sig,czero,nmax,plan,nn,FFTW_FORWARD,FFTW_ESTIMATE,dom,dt,nk,npts,fre,seisout,comb
implicit none
integer nerr,i
real(8) amax,f
real(8) ss1(nmax),ss2(nmax),ss3(nmax)
complex(8) s1(nmax),s2(nmax),s3(nmax)
complex(8) sa(nmax),sb(nmax),sc(nmax)
complex(8) saf(nmax),sbf(nmax),scf(nmax)
!
sa=czero
sb=czero
s1=czero
s2=czero
saf=czero
sbf=czero
seisout=czero
sa(1:npts)=dcmplx(dble(sig(1:npts,1)),0.d0)
sb(1:npts)=dcmplx(dble(sig(1:npts,2)),0.d0)
call dfftw_plan_dft_1d(plan,nn,sa,saf,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_1d(plan,nn,sb,sbf,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
call smoothf(2,saf,ss1)
call smoothf(2,sbf,ss2) ! smooth the spectrum
do i=1,nk 
   amax=max(ss1(i),ss2(i))
   f = dble((i-1))*dom
   if( f .ge. dble(fre(1,2)) .and. f .le. dble(fre(4,2)) ) then
      if(comb.eq.1)then
         s1(i)=saf(i)/amax
         s2(i)=sbf(i)/amax
      else
         s1(i)=saf(i)/ss1(i)
         s2(i)=sbf(i)/ss2(i)
      endif
   endif
enddo
! taper the spectrum
call  taperf(s1,2)
call  taperf(s2,2)
seisout(1:nn,1) = s1(1:nn)
seisout(1:nn,2) = s2(1:nn)
end subroutine
