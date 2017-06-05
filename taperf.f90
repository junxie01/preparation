! Tapering subroutine itself in frequency domain
! 2017/06/04
subroutine taperf(sf,id)
use globe_data
integer   i,j,id
real(8)   f,d1,d2,dpi,ss,s(nmax)
complex(8)   sf(nmax)
! ---
dpi = datan(1.0d0)*4.0d0
s=0d0
do i = 1,nk
   f = dble(i-1)*dom
   if(f.le.dble(fre(1,id))) then
      cycle
   else if(f.le.dble(fre(2,id))) then
      d1 = dpi/(dble(fre(2,id))-dble(fre(1,id)))
      ss = 1.0d0
      do j = 1,npow
         ss = ss*(1.d0-dcos(d1*(dble(fre(1,id))-f)))/2.0d0
      enddo
      s(i) = ss
   else if(f.le.dble(fre(3,id))) then
      s(i) = 1.0d0
   else if(f.le.dble(fre(4,id))) then
      d2 = dpi/(dble(fre(4,id))-dble(fre(3,id)))
      ss = 1.0d0
      do j = 1,npow
         ss = ss*(1.0d0+dcos(d2*(fre(3,id)-f)))/2.0d0
      enddo
      s(i) = ss
   endif
enddo
do i=1,nk
   sf(i) = sf(i)*s(i)
enddo
return
end subroutine
