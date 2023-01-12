
subroutine rtmx_3 (filnam,wt_par,f41)

!*************************************************************************
!* RTMX function in MTG2 version 2.21 
!* Estimation of kernel matrices
!* Author: S. Hong Lee, University of South Australia
!* (C) 2021-2022 NIAS, RDA - All Rights Reserved under GNU LGPL
!*************************************************************************
 
implicit none
real::sclf(100)     ! scale factor, check if pops # > 100
integer::nm,nid,i,j,k,zi,zj,im,xi,yi,wt,npop,eoro
integer,allocatable::k1(:)

character(len=100),allocatable::pedor(:,:)
character(len=50),allocatable::markt(:)
real,allocatable::mafreq2(:,:),mpos(:),htz(:,:),htz1(:)
!real,allocatable::kvm(:,:),kvm2(:,:)
double precision,allocatable::kvm(:,:),kvm2(:,:),rtmx(:,:),rtmx2(:,:),rtmx_off(:,:)
real,allocatable::rtmx_tot(:,:)

character(len=1),allocatable::genom(:,:),mafreq1(:,:)
character(len=3)::cdum4
character(len=3),allocatable::cno(:)

real::x1,x2,x3,htz2a,htz2b
integer,allocatable::wt_idx1(:,:),nm_idx(:)
double precision::diags,v1,v2
double precision,allocatable::diag1(:)

!sequential storing
integer*1,allocatable::sqgenom(:,:),sex(:) 
integer*1 b1,igen(0:1,0:1),w1,w2

logical::wt_logic

!progressing report **************************
character*3 bsp

!reading command line *****************
integer::narg,io
character(len=64)::filnam,cdum,f41,wt_par,wt_file

!real time **************************
character(8)::date
character(10)::time
character(5)::zone

call date_and_time(date,time,zone)
print'(a,2x,a,2x,a,2x,a)',"== input",date,time,zone


nid=0
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
do
  nid=nid+1
  read(35,*,iostat=io)cdum4
  if (io.ne.0) exit
end do
close(35)
nid=nid-1
print*,'no. ID    : ',nid

nm=0
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file
do
  nm=nm+1
  read(36,*,iostat=io)cdum4
  if (io.ne.0) exit
end do
close(36)
nm=nm-1
print*,'no. marker: ',nm


ALLOCATE(mpos(nm),cno(nm),mafreq1(nm,2)) !,mafreq2(nm,2)) 
allocate(pedor(nid,4))
allocate(sqgenom(int(nid/4)+1,nm),sex(nid))
!allocate(htz(nm,2),markt(nm))
allocate(markt(nm),nm_idx(nm))


if (wt_par.ne."null") then
  open (unit=51,file=wt_par,status='old')
  read(51,*)npop
  read(51,*)wt_file
  do i=1,npop
    read(51,*)sclf(i)
  end do
  read(51,*)eoro       !1=2pq or 2=var(x)
  close(51)
else
  print*,'pop info file, scale factors for pops. >>> check'
  pause 
end if
!print*,npop,sclf,eoro

allocate(wt_idx1(npop,nid),k1(npop),diag1(npop))   !,wt_idx2(nid))
allocate(mafreq2(nm,npop),htz(nm,npop),htz1(npop))
inquire(file=wt_file,exist=wt_logic)
if (wt_logic) then
  open (unit=39,file=wt_file,status='old')      !pop info file
  k1=0 !;k2=0
  do i=1,nid
    read(39,*)cdum4,cdum4,wt             !pop info, 1, 2,...,
    k1(wt)=k1(wt)+1
    wt_idx1(wt,k1(wt))=i
    !print*,wt,k1(wt),wt_idx1(wt,k1(wt))
  end do
  print*,'pop info file:',wt_file
  print*,''
else
  print*,'file missing >>>',wt_file
  pause
end if
do i=1,npop
  print*,'sample size for pop:',i,k1(i)
end do
print*,''



open (unit=38,file=trim(filnam)//'.bed',status='old',access='stream',form='unformatted')
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file


do i=1,nm
  read(36,*)cno(i),markt(i),cdum4,mpos(i),mafreq1(i,1),mafreq1(i,2)
end do


!progressing report***************************************************
bsp(1:1)=char(8); bsp(2:2)=char(8); bsp(3:3)=char(8)

!for fam file (ped) ***************************************************
do i=1,nid
  read(35,*)pedor(i,:),sex(i),cdum4
end do !***************************************************************


! for pedigree and genotype (PLINK .bed, * files) ****************************
read(38)b1           !plink magic number 1


if (.not.btest(b1,0).and..not.btest(b1,1).and.btest(b1,2).and.btest(b1,3).and.&
& .not.btest(b1,4).and.btest(b1,5).and.btest(b1,6).and..not.btest(b1,7)) then
  write(*,'(a7)',advance='no')' 1 - ok'
else
  print*,'this may not be PLINK .bed file - check >>>'
end if

read(38)b1           !plink magic number 2
if (btest(b1,0).and.btest(b1,1).and..not.btest(b1,2).and.btest(b1,3).and.&
& btest(b1,4).and..not.btest(b1,5).and..not.btest(b1,6).and..not.btest(b1,7)) then
  write(*,'(a7)',advance='no')' 2 - ok'
else
  print*,'this may not be PLINK .bed file - check >>>'
end if

read(38)b1           !mode 

if (btest(b1,0)) then
  write(*,'(a37)')'  SNP-major mode for PLINK .bed file'
else
  print*,'should not be individual mode - check >>>'
end if


do im=1,nm
  if (nm<100 .or. mod(im,int(nm/100))==0) then
    write (*,'(a3,i3)',advance='no')bsp,int((.1*im)/(.1*nm)*100)
    !print*,im
!20  format (A3,I8)
  end if

  do i=1,int((nid-1)/4)+1    !no. ID / 4
    read(38)sqgenom(i,im)
  end do
  !print*,igenom(i,1:20)
  !pause

end do
write (*,'(a6)')'% done'    !closing progress report

close(35)
close(38)
close(36)

!print*,k1,k2

!check duplications in ped file **********************************
do i=1,nid
  do j=i+1,nid
    if (pedor(i,1)==pedor(j,1).and.pedor(i,2)==pedor(j,2)) then
      print*,'ID duplications in .ped: ',trim(pedor(i,1)),trim(pedor(i,2))
      pause
    end if
  end do
end do !**********************************************************

call date_and_time(date,time,zone)
print'(a,2x,a,2x,a,2x,a)',"== start",date,time,zone

!indexed genotype coefficients into igen********************************
igen(0,0)=2
igen(0,1)=1
igen(1,1)=0
igen(1,0)=3           !missing *****************************************

PRINT*,''

do zi=1,npop
  PRINT*,'allele frequency for pop',zi
  !$OMP PARALLEL DO PRIVATE (yi,xi,w1) 
  do im=1,nm
    x1=0;x2=0;x3=0
    do j=1,k1(zi)  !nid
      i=wt_idx1(zi,j)
      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom(yi,im),xi*2,1),ibits(sqgenom(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        x1=x1+w1
        x2=x2+1
        x3=x3+w1**2
      end if
    end do
    mafreq2(im,zi)=x1/x2     !mean*2
    if (eoro==2) then
      htz(im,zi)=(x3-(x1**2/x2))/(x2)
    elseif (eoro==1) then
      htz(im,zi)=mafreq2(im,zi)*(1-mafreq2(im,zi)/2)
    end if
    !print*,zi,im,mafreq2(im,zi),htz(im,zi)
  end do
  !$OMP END PARALLEL DO
end do

allocate(rtmx_tot(nid,nid))

do zi=1,npop
  allocate(kvm(k1(zi),nm))
  allocate(rtmx(k1(zi),k1(zi)))

  print*,'rtmx estimation - scale factor for pop:',zi,sclf(zi)
  !$OMP PARALLEL DO PRIVATE (yi,xi,w1,x3) 
  do im=1,nm
    do j=1,k1(zi)  !nid
      i=wt_idx1(zi,j)
      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom(yi,im),xi*2,1),ibits(sqgenom(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        kvm(j,im)=(w1-mafreq2(im,zi))*(htz(im,zi)**(0.5*sclf(zi)))
        !print*,zi,im,j,kvm(i,im)
      else
        kvm(j,im)=0
      end if
      !if (j==1) print'(10000i1)',w1
    end do
    !if (im==1) print*,mafreq2(im,1),htz(im,1),0.5*sclf(1)
    !if (im==1) print*,mafreq2(im,2),htz(im,2),0.5*sclf(2)
  end do
  !$OMP END PARALLEL DO

  htz1(zi)=0
  do im=1,nm
    if (htz(im,zi) > 0 .and. htz(im,zi) < 1) then  !check, exclude NaN htz
      htz1(zi)=htz1(zi)+(htz(im,zi)**(1+sclf(zi)))/k1(zi)   ! for correction
    end if
    !print*,htz1(zi)
  end do
  !pause

  call dgemm ('N','T',k1(zi),k1(zi),nm,1.0D0,kvm,k1(zi),kvm,k1(zi),0.0D0,rtmx,k1(zi))
  ! correction (see eq. 9 in Momin's paper)
  rtmx=rtmx+htz1(zi)
  diag1(zi)=0
  do i=1,k1(zi)
    diag1(zi)=diag1(zi)+rtmx(i,i)
  end do
  diag1(zi)=diag1(zi)/k1(zi)

  do zj=1,(zi-1)
    !print*,zj,zi
    allocate(kvm2(k1(zj),nm))
    allocate(rtmx_off(k1(zi),k1(zj)))

    !print*,'rtmx estimation - scale factor for pop:',zj,sclf(zj)
    !$OMP PARALLEL DO PRIVATE (yi,xi,w1) 
    do im=1,nm
      do j=1,k1(zj)  !nid
        i=wt_idx1(zj,j)
        yi=int((i-1)/4)+1
        xi=mod((i-1),4)
        w1=igen(ibits(sqgenom(yi,im),xi*2,1),ibits(sqgenom(yi,im),xi*2+1,1))
        if (w1.ne.3) then
          kvm2(j,im)=(w1-mafreq2(im,zj))*(htz(im,zj)**(0.5*sclf(zj)))
        else
          kvm2(j,im)=0
        end if
        !if (j==1) print'(10000i1)',w1
      end do
    end do
    !$OMP END PARALLEL DO

    call dgemm ('N','T',k1(zi),k1(zj),nm,1.0D0,kvm,k1(zi),kvm2,k1(zj),0.0D0,rtmx_off,k1(zi))

    ! correction (see eq. 9 in Momin's paper)
    k=0;htz2a=0;htz2b=0
    do im=1,nm
      if (htz(im,zi)>0.and.htz(im,zi)<1.and.htz(im,zj)>0.and.htz(im,zj)<1) then  !check, exclude NaN htz
        !htz2=htz2+sqrt((htz(im,zi)**(1+sclf(zi)))/k1(zi)) * ((htz(im,zj)**(1+sclf(zj)))/k1(zj))   ! for correction
        htz2a=htz2a+((htz(im,zi)**(1+sclf(zi)))/k1(zi)) 
        htz2b=htz2b+((htz(im,zj)**(1+sclf(zj)))/k1(zj))   ! for correction
        k=k+1
        nm_idx(k)=im
      end if
      !print*,htz2a,htz2b
    end do
    
    rtmx_off=rtmx_off+sqrt(htz2a*htz2b)
    !print*,rtmx_off

    !between pops
    v1=0;v2=0
    do i=1,k1(zi)
      do im=1,k
        v1=v1+kvm(i,nm_idx(im))**2
      end do
      v1=v1+htz2a
    end do
    do j=1,k1(zj)
      do im=1,k
        v2=v2+kvm2(j,nm_idx(im))**2
      end do
      v2=v2+htz2b
    end do
    !diags=diag1(zi)*diag1(zj)
    diags=(v1/real(k1(zi)))*(v2/real(k1(zj)))
    !print*,v1/real(k1(zi)),v2/real(k1(zj))

    !$OMP PARALLEL DO PRIVATE (i,j) 
    do i=1,k1(zi)
      do j=1,k1(zj)
        !rtmx_tot(wt_idx1(zi,i),wt_idx1(zj,j))=real(rtmx_off(i,j))
        !rtmx_tot(wt_idx1(zj,j),wt_idx1(zi,i))=real(rtmx_off(i,j))
        rtmx_tot(wt_idx1(zi,i),wt_idx1(zj,j))=real(rtmx_off(i,j)/sqrt(diags))
        rtmx_tot(wt_idx1(zj,j),wt_idx1(zi,i))=real(rtmx_off(i,j)/sqrt(diags))
      end do
    end do
    !$OMP END PARALLEL DO
    deallocate(kvm2)
    deallocate(rtmx_off)

  end do !zj

  !$OMP PARALLEL DO PRIVATE (i,j) 
  do i=1,k1(zi)
    do j=1,i
      rtmx_tot(wt_idx1(zi,i),wt_idx1(zi,j))=real(rtmx(i,j)/diag1(zi))
      rtmx_tot(wt_idx1(zi,j),wt_idx1(zi,i))=real(rtmx(i,j)/diag1(zi))
    end do
  end do
  !$OMP END PARALLEL DO

  deallocate(kvm)
  deallocate(rtmx)

end do  !zi

!print*,diag1
!print*,diag2
!print*,sqrt(diag1*diag2)
!print*,rtmx_off(1,1:5)
!print*,rtmx_off(2,1:5)
!print*,rtmx_off(3,1:5)
!print*,rtmx_off(4,1:5)
!print*,rtmx_off(5,1:5)
!print*,''

!print'(10000f14.8)',kvm(1,:)
!print'(10000f14.8)',kvm2(1,:)

call date_and_time(date,time,zone)
print'(a,2x,a,2x,a,2x,a)',"== finish",date,time,zone
print*,''

open (UNIT=41,FILE=trim(f41),STATUS='unknown')
do i=1,nid
  do j=1,i 
    write(41,'(i8,i8,f14.8)')i,j,real(rtmx_tot(i,j))
  end do
end do
close(41)
print*,'Relationship matrix (RTMX) stored in  ',trim(f41)

call date_and_time(date,time,zone)
print'(a,2x,a,2x,a,2x,a)',"== output",date,time,zone

end subroutine

