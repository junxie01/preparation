program pre_processing 
! do 3 component preprocessing 
! hope it can process one or two component too
! cut the data to avoid the side problems due to bp filter
! 2017/06/04 --introduce globe_data
use sacio
use globe_data
implicit none
type(sac_head):: sachead(3)
logical ext
real stla,stlo,evla,evlo,t1
character (2)net(nstmax)
character (3)com(3),co
character (20)year_day,nd
character (10)sta(nstmax),kst
character (180)sacf(3),command
character (180)out_w(3),out(3)
character (180)dir_data,inpar,sta_list
integer nst,ist,itest
integer iy,id,is
integer nerr,iseg
integer jday,icheck,n1
integer begday,endday,ic
integer multpt,dsec,dseg,nseg
integer nzhour,nzmin,nzsec,nzmsec
integer year_b,year_e,day_b,day_e

if (iargc().ne.1)then
   write(*,*)'Usage: preparation param.dat '
   write(*,*)'param.dat:'
   write(*,*)'staion_list'
   write(*,*)'year_b day_b year_e day_e'
   write(*,*)'dsec multpt com ncom comb'
   write(*,*)'fa fb f1 f2'
   write(*,*)'t1, npts'
   write(*,*)'data_dir'
   call exit(-1)
endif
call getarg(1,inpar)
inquire(file=inpar,exist=ext)
if(.not.ext)stop "Hey dude, forget something?"
! read the control parameter
open(10,file=inpar)
read(10,'(a180)')sta_list
read(10,*)year_b,day_b,year_e,day_e
read(10,*)dsec,multpt,co,ncom,comb
read(10,*)fre(2,1),fre(3,1),fre(2,2),fre(3,2)
read(10,*)t1,npts
read(10,'(a180)')dir_data
close(10)

! define the component Z,N,E
com(1)=trim(co)//'Z'
com(2)=trim(co)//'N'
com(3)=trim(co)//'E'
if(ncom.eq.1)com(1)=trim(co)
if(ncom.eq.2)then
   com(1)=trim(co)//'N'
   com(2)=trim(co)//'E'
endif
! define the corner frequency
! you can define them by yourself
fre(1,1)=0.85*fre(2,1)
fre(4,1)=1.15*fre(3,1)
fre(1,2)=0.85*fre(2,2)
fre(4,2)=1.15*fre(3,2)

write(*,*)1.0/fre(4,1),1.0/fre(3,1),1.0/fre(2,1),1.0/fre(1,1)
write(*,*)1.0/fre(4,2),1.0/fre(3,2),1.0/fre(2,2),1.0/fre(1,2)
nn=2
npow=1
do while(nn.lt.npts)
   nn=nn*2
   npow=npow+1
enddo
nk = nn/2+1

! read in station list
open(11,file=sta_list)          
do ist=1,nstmax
   read(11,*,err=13,end=13) net(ist),sta(ist)
enddo
13 close(11)
nst=ist-1

dseg=int((1-real(multpt)/100.0)*dsec) ! the left points without overlapping
nseg=int((86400-dsec)/dseg)+1         ! number of segments per day.
do iy=year_b,year_e                   ! loop over year
   jday=365
   if(mod(iy,4).eq.0.and.mod(iy,100).ne.0.or.mod(iy,400).eq.0)jday=366
   endday=day_e
   if(iy.ne.year_e)endday=jday
   begday=day_b
   if(iy.ne.year_b)begday=1
   do id=begday,endday                ! loop over day
      write(year_day,'(i0,"_",i3.3)')iy,id
      write(*,'("Deal with ",i0,"/",i3.3)')iy,id
      do is=1,nst                     ! loop over station 
         do iseg=1,nseg               ! loop over segments
            nzhour=(iseg-1)*dseg/3600
            nzmin=mod((iseg-1)*dseg,3600)/60
            nzsec=mod(mod((iseg-1)*dseg,3600),60)
            icheck=0
            do ic=1,ncom              ! loop over components
               write(sacf(ic),'(1a,"/",1a,"/",1a,"_",i2.2,"_",i2.2,"_",i2.2,"_",1a,"_",1a,"_",1a,".SAC")')&
               trim(dir_data),trim(year_day),trim(year_day),nzhour,nzmin,nzsec,trim(net(is)),trim(sta(is)),trim(com(ic))
               inquire(file=sacf(ic),exist=ext)
               icheck=ic
            enddo
            sig=0
            if(icheck.ne.ncom)cycle                   ! check whether all components exist
            do ic=1,ncom                              ! read all data
               call read_sachead(trim(sacf(ic)),sachead(ic),nerr)
               if(ic.eq.1)n1=int(t1/sachead(ic)%delta)+1
               call read_sac(trim(sacf(ic)),sigt(:),sachead(ic),nerr)
               if(nerr.eq.-1)exit
               sig(1:npts,ic)=sigt(n1:n1+npts-1)
            enddo
            !do itest=1,npts
            !   write(*,*)sig(itest,1)
            !enddo
            if(nerr.eq.-1)cycle                       ! read sac file fails
            dt=sachead(1)%delta
            dom=dble(1.0/nn/dt)
            nf=int(0.005d0/dom)
            halfn=int(64/dt) 
            !write(*,*)'Do time domain normalisation'
            if(ncom.eq.1)call do_norm1
            if(ncom.eq.2)call do_norm2
            if(ncom.eq.3)call do_norm3
            !write(*,*)'Do frequency domain whitening'
            if(ncom.eq.1)call do_whiten1
            if(ncom.eq.2)call do_whiten2
            if(ncom.eq.3)call do_whiten3
            do ic=1,ncom
               sachead(ic)%npts=nn
               !sachead(ic)%delta=sngl(dom)
               sigo(1:nn,ic)=sngl(dreal(seisout(1:nn,ic)))
               call write_sac(trim(sacf(ic))//".re",sigo(:,ic),sachead(ic),nerr)
               sigo(1:nn,ic)=sngl(dimag(seisout(1:nn,ic)))
               call write_sac(trim(sacf(ic))//'.im',sigo(:,ic),sachead(ic),nerr)
               !sigo(1:nn,ic)=sngl(dsqrt(dimag(seisout(1:nn,ic))**2+dreal(seisout(1:nn,ic))**2))
               !call write_sac(trim(sacf(ic))//'.amp',sigo(:,ic),sachead(ic),nerr)
            enddo
            !write(*,*)'Write done!'
         enddo      ! end loop over segments
      enddo         ! end loop over station
   enddo            ! end loop over day
enddo               ! end loop over year
write(*,*)"All done for one station pre-processing!"
end program
