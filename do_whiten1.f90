!2017/06/04
subroutine do_whiten1
        use globe_data,only: sig,czero,nmax,plan,nn,FFTW_FORWARD,FFTW_ESTIMATE,dom,dt,nk,npts,fre,seisout
implicit none
integer i
real(8) amax,f
real(8) ss1(nmax)
complex(8) s1(nmax)
complex(8) sa(nmax)
complex(8) saf(nmax)
sa=czero
s1=czero
saf=czero
sa(1:npts)=dcmplx(dble(sig(1:npts,1)),0.d0)
call dfftw_plan_dft_1d(plan,nn,sa,saf,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
saf(1) = saf(1)/2.0d0
saf(nk) = dcmplx(dreal(saf(npts)),0.0d0)
saf(1:nk)=saf(1:nk)*dble(dt)

call smoothf(2,saf,ss1) ! smooth the spectrum
do i=1,nk 
   amax=ss1(i)
   f = dble((i-1))*dom
   if( f .ge. dble(fre(1,2)) .and. f .le. dble(fre(4,2)) ) then
      s1(i)=saf(i)/amax
   endif
enddo
call taperf(s1,2)
seisout(1:nn,1) = s1(1:nn)
end subroutine
