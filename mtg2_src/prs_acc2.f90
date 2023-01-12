
subroutine prs_acc (filnam,gwasf,targ,fl4,f41)

!*************************************************************************
!* var(beta), var(PRS), var(R2) function in MTG2 version 2.22 
!* Author: S. Hong Lee, University of South Australia
!* (C) 2021-2022 NIAS, RDA - All Rights Reserved under GNU LGPL
!*************************************************************************
 
implicit none
real::sclf(100)     ! scale factor, check if pops # > 100
integer::nm,nm2,nid,i,j,k,zi,zj,im,xi,yi,wi,npop,eoro,intv,totn,mn,tn
integer,allocatable::k1(:),startv(:),endv(:)

character(len=100),allocatable::pedor(:,:),tpedor(:,:)
character(len=50),allocatable::markt(:),tmarkt(:)
real,allocatable::mafreq2(:,:),mpos(:),htz(:,:)
real,allocatable::tmafreq2(:,:),tmpos(:),thtz(:,:)
double precision,allocatable::kvm(:,:),kvm2(:,:),rtmx(:,:),snpm11(:,:)
double precision,allocatable::tkvm(:,:),tkvm2(:,:)

character(len=1),allocatable::genom(:,:),mafreq1(:,:),tmafreq1(:,:)
character(len=3)::cdum4,cdum
character(len=3),allocatable::cno(:),tcno(:)

real::x1,x2,x3,h2
double precision::LA,signs
double precision,allocatable::beta(:),beta_blup(:),beta_blup_v(:,:)
double precision,allocatable::beta2(:),beta_blup2(:),beta_blup_v2(:,:)
double precision,allocatable::blup(:),blup_v(:,:)

!sequential storing
integer*1,allocatable::sqgenom(:,:),sex(:) 
integer*1,allocatable::tsqgenom(:,:),tsex(:) 
integer*1 b1,igen(0:1,0:1),w1,w2

!progressing report **************************
character*3 bsp

!reading command line *****************
integer::narg,io
character(len=64)::filnam,f41,gwasf,fl4,targ

double precision::vcva(1,1,1),vcve(1,1)
logical::file_exist

  tn=1;mn=1     !check
  inquire(file=fl4,exist=file_exist)
  if (file_exist) then
    open (UNIT=44,FILE=fl4,STATUS='old')
    do i=1,tn
      read(44,*,iostat=io)cdum,vcve(i,i)
      if (io.ne.0) exit
    end do
    do wi=1,mn
      do xi=1,tn
        read(44,*,iostat=io)cdum,vcva(wi,xi,xi)
        if (io.ne.0) exit
      end do
      do xi=1,tn
        do yi=1,xi-1
          read(44,*,iostat=io)cdum,vcva(wi,xi,yi)
          if (io.ne.0) exit
        end do
      end do
    end do
    close(44)
  else
    print*,'starting value have to be provided for this BLUP function >>>'
    pause
  end if

  h2=vcva(1,1,1)/(vcva(1,1,1)+vcve(1,1))

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

!set a signle pop, scale factor=-1 (alpha=-0.5), uisng var(x)
npop=1; sclf=-1; eoro=2

allocate(mpos(nm),cno(nm),mafreq1(nm,2)) !,mafreq2(nm,2)) 
allocate(pedor(nid,4))
allocate(sqgenom(int(nid/4)+1,nm),sex(nid))
allocate(markt(nm))
allocate(mafreq2(nm,npop),htz(nm,npop))

!by chromosomes and by blocks
intv=100 !nm !100
if (mod(nm,intv)==0) then
  totn=nm/intv
else
  totn=(nm/intv)+1
end if
allocate(startv(totn),endv(totn))

do i=1,totn
  startv(i)=i*intv-intv+1
  endv(i)=i*intv
  print*,startv(i),endv(i)
end do
endv(totn)=nm


open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file
do i=1,nm
  read(36,*)cno(i),markt(i),cdum4,mpos(i),mafreq1(i,1),mafreq1(i,2)
end do
close(36)

open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
!for fam file (ped) ***************************************************
do i=1,nid
  read(35,*)pedor(i,:),sex(i),cdum4
end do !***************************************************************
close(35)

!progressing report***************************************************
bsp(1:1)=char(8); bsp(2:2)=char(8); bsp(3:3)=char(8)

open (unit=38,file=trim(filnam)//'.bed',status='old',access='stream',form='unformatted')
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
close(38)


!check duplications in ped file **********************************
do i=1,nid
  do j=i+1,nid
    if (pedor(i,1)==pedor(j,1).and.pedor(i,2)==pedor(j,2)) then
      print*,'ID duplications in .ped: ',trim(pedor(i,1)),trim(pedor(i,2))
      pause
    end if
  end do
end do !**********************************************************

!indexed genotype coefficients into igen********************************
igen(0,0)=2
igen(0,1)=1
igen(1,1)=0
igen(1,0)=3           !missing *****************************************

PRINT*,''

allocate(beta(nm),beta_blup(nm),beta_blup_v(nm,nm))
!allocate(blup(nid),blup_v(nid,nid),kvm(nid,nm))
allocate(kvm(nid,nm))
beta_blup_v=0

!reading GWAS summary
!open (unit=42,file="plink.qassoc",status='old')      !fam file
!open (unit=42,file="toy.qassoc",status='old')      !fam file
!open (unit=42,file="toy2.qassoc",status='old')      !fam file
!open (unit=42,file="discovery2.qassoc",status='old')      !fam file
open (unit=42,file=trim(gwasf),status='old')      !fam file
read(42,*)cdum4
do i=1,nm
  read(42,*)cdum4,cdum4,cdum4,cdum4,signs,cdum4,beta(i)
  beta(i)=sqrt(beta(i))
  if (signs < 0.0D0) beta(i)=beta(i)*-1.0D0
  !print*,beta(i)
end do
close(42)

  !PRINT*,'allele frequency for pop',zi
  !$OMP PARALLEL DO PRIVATE (yi,xi,w1) 
  do im=1,nm
    x1=0;x2=0;x3=0
    do j=1,nid
      i=j !wt_idx1(zi,j)
      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom(yi,im),xi*2,1),ibits(sqgenom(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        x1=x1+w1
        x2=x2+1
        x3=x3+w1**2
      end if
    end do
    mafreq2(im,1)=x1/x2     !mean*2
    if (eoro==2) then
      htz(im,1)=(x3-(x1**2/x2))/(x2)
    elseif (eoro==1) then
      htz(im,1)=mafreq2(im,1)*(1-mafreq2(im,1)/2)
    end if
  end do
  !$OMP END PARALLEL DO


do zi=1,totn
  nm2=endv(zi)-startv(zi)+1
  !allocate(kvm(nid,nm),kvm2(nid,nm),rtmx(nm,nm))
  allocate(kvm2(nid,nm2),rtmx(nm2,nm2),snpm11(nm2,nm2))
  allocate(beta_blup2(nm2),beta_blup_v2(nm2,nm2))

  !$OMP PARALLEL DO PRIVATE (yi,xi,w1,x3) 
  !do im=1,nm
  do im=startv(zi),endv(zi)
    do i=1,nid
      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom(yi,im),xi*2,1),ibits(sqgenom(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        !kvm(i,im)=(w1-mafreq2(im,1))*(htz(im,1)**(0.5*sclf(1)))
        kvm2(i,(im-startv(zi)+1))=(w1-mafreq2(im,1))*(htz(im,1)**(0.5*sclf(1)))
      else
        !kvm(i,im)=0
        kvm2(i,(im-startv(zi)+1))=0
      end if
      !print*,kvm(i,im)
    end do
    !pause
  end do
  !$OMP END PARALLEL DO
  !call dgemm ('T','N',nm,nm,nid,1.0D0,kvm,nid,kvm,nid,0.0D0,rtmx,nm)
  call dgemm ('T','N',nm2,nm2,nid,1.0D0,kvm2,nid,kvm2,nid,0.0D0,rtmx,nm2)

  !h2=0.5
  snpm11=rtmx*(h2/(nid*nm2))
  do i=1,nm2
    snpm11(i,i)=snpm11(i,i)+(1-h2)/nid
  end do

  call cholesky_inv (snpm11,nm2,LA)
  print*,nm2,LA

  !SNP BLUP
  !call dgemm ('N','N',nm2,1,nm2,dble(h2/nm),snpm11,nm2,beta,nm,0.0D0,beta_blup,nm)
  call dgemm ('N','N',nm2,1,nm2,dble(h2/nm),snpm11,nm2,beta(startv(zi):endv(zi)),nm2,0.0D0,beta_blup2,nm2)
  !var(SNP BLUP)
  !call dgemm ('N','N',nm,nm,nm,dble(h2/nm)**2,snpm11,nm,rtmx/nid,nm,0.0D0,beta_blup_v,nm)
  call dgemm ('N','N',nm2,nm2,nm2,dble(h2/nm)**2,snpm11,nm2,rtmx/nid,nm2,0.0D0,beta_blup_v2,nm2)


  beta_blup(startv(zi):endv(zi))=beta_blup2
  !do i=1,nm2
  !  beta_blup_v(startv(zi)-1+i)=beta_blup_v2(i,i)
  !end do
  beta_blup_v(startv(zi):endv(zi),startv(zi):endv(zi))=beta_blup_v2
  kvm(:,startv(zi):endv(zi))=kvm2

  deallocate(kvm2,rtmx,snpm11,beta_blup2,beta_blup_v2)
end do !zi


print*,"target data"
!filnam="target2"
filnam=trim(targ)
!nm=2872 !2877   !check
!nid=200 !1000

nid=0
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
do
  nid=nid+1
  read(35,*,iostat=io)cdum4
  if (io.ne.0) exit
end do
close(35)
nid=nid-1
print*,'no. ID in target set: ',nid

nm=0
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file
do
  nm=nm+1
  read(36,*,iostat=io)cdum4
  if (io.ne.0) exit
end do
close(36)
nm=nm-1
print*,'no. marker in target set: ',nm

allocate(tmpos(nm),tcno(nm),tmafreq1(nm,2)) !,mafreq2(nm,2)) 
allocate(tpedor(nid,4))
allocate(tsqgenom(int(nid/4)+1,nm),tsex(nid))
allocate(tmarkt(nm))
allocate(tmafreq2(nm,npop),thtz(nm,npop))

open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file
do i=1,nm
  read(36,*)tcno(i),tmarkt(i),cdum4,tmpos(i),tmafreq1(i,1),tmafreq1(i,2)
end do
close(36)

open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
!for fam file (ped) ***************************************************
do i=1,nid
  read(35,*)tpedor(i,:),tsex(i),cdum4
end do !***************************************************************
close(35)

open (unit=38,file=trim(filnam)//'.bed',status='old',access='stream',form='unformatted')
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
  do i=1,int((nid-1)/4)+1    !no. ID / 4
    read(38)tsqgenom(i,im)
  end do
  !print*,igenom(i,1:20)
  !pause
end do
close(38)

allocate(tkvm2(nid,nm),tkvm(nid,nm),blup(nid),blup_v(nid,nid))

  !PRINT*,'allele frequency for pop',zi
  !$OMP PARALLEL DO PRIVATE (yi,xi,w1) 
  do im=1,nm
    x1=0;x2=0;x3=0
    do j=1,nid
      i=j !wt_idx1(zi,j)
      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(tsqgenom(yi,im),xi*2,1),ibits(tsqgenom(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        x1=x1+w1
        x2=x2+1
        x3=x3+w1**2
      end if
    end do
    tmafreq2(im,1)=x1/x2     !mean*2
    if (eoro==2) then
      thtz(im,1)=(x3-(x1**2/x2))/(x2)
    elseif (eoro==1) then
      thtz(im,1)=tmafreq2(im,1)*(1-tmafreq2(im,1)/2)
    end if
  end do
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO PRIVATE (yi,xi,w1,x3) 
  do im=1,nm
    do i=1,nid
      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(tsqgenom(yi,im),xi*2,1),ibits(tsqgenom(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        tkvm(i,im)=(w1-tmafreq2(im,1))*(thtz(im,1)**(0.5*sclf(1)))
      else
        tkvm(i,im)=0
      end if
    end do
  end do
  !$OMP END PARALLEL DO

print*,"EBV"
call dgemm ('N','N',nid,1,nm,1.0D0,tkvm,nid,beta_blup,nm,0.0D0,blup,nid)
print*,"var(EBV)"
call dgemm ('N','N',nid,nm,nm,1.0D0,tkvm,nid,beta_blup_v,nm,0.0D0,tkvm2,nid)
call dgemm ('N','T',nid,nid,nm,1.0D0,tkvm2,nid,tkvm,nid,0.0D0,blup_v,nid)


open (UNIT=41,FILE=trim(f41)//".snpblup",STATUS='unknown')
do i=1,nm
!  do j=1,i 
!    write(41,'(i8,i8,f14.8)')i,j,real(rtmx(i,j)/nid)
!    !write(41,'(i8,i8,f14.8)')i,j,real(rtmx(i,j))
!  end do
  write(41,'(i8,i8,f14.8,f14.8)')i,i,real(beta_blup(i)),real(sqrt(beta_blup_v(i,i)))
end do
close(41)

open (UNIT=41,FILE=trim(f41)//".indblup",STATUS='unknown')
do i=1,nid
  write(41,'(i8,i8,f14.8,f14.8)')i,i,real(blup(i)),real(blup_v(i,i))**0.5
end do
close(41)
!print*,'SNP correlation matrix stored in  ',trim(f41)

end subroutine

