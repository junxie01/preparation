! 2017/06/04
module globe_data
integer,parameter :: nmax=4000000,nstmax=1000
real(4),dimension(4,2):: fre
real(4),dimension(nmax,3):: sig,sigo
real(4),dimension(nmax) :: sigt
real :: dt
real(8) :: dom
integer :: halfn
integer :: ncom,comb,npts
integer :: nn,npow,nk,nf
complex(8),dimension(nmax,3):: seisout
complex(8),parameter:: czero=(0.0d0,0.0d0)
!integer,parameter :: FFTW_ESTIMATE=0,FFTW_MEASURE=1
integer,parameter :: FFTW_ESTIMATE=64,FFTW_MEASURE=1
integer,parameter :: FFTW_FORWARD=-1,FFTW_BACKWARD=1
integer(8) :: plan,plan1,plan2,plan3

end module
