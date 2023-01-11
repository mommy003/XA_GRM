
subroutine h2_mtpreg (filnam,yv,nid,fl5)

!Estimating SNP-h2 using multiple regression 
!Sang Hong Lee (Apr/22)
 
implicit none
integer::nm,nid,sn,dn,i,j,k,k2,zi,zj,im,ma1,idum,nm2,nm3,xi,yi,xi2,yi2

character(len=100),allocatable::pedor(:,:)
character(len=50),allocatable::markt(:)
real,allocatable::mafreq2(:),mpos(:),htz(:)
double precision,allocatable::kvm(:,:),xx(:,:),xy(:,:),sol(:,:),yhat(:,:)
double precision::LA,tyv(nid,1)

character(len=1),allocatable::genom(:,:),mafreq1(:,:)
character(len=3)::cdum4
character(len=3),allocatable::cno(:)

real::x1,x2,x3,yv(nid,1),sse,sst

!sequential storing
integer*1,allocatable::sqgenom(:,:),sex(:) !,igenom(:,:)
integer*1 b1,igen(0:1,0:1),w1,w2

!progressing report **************************
character*3 bsp

!reading command line *****************
integer::narg,io
character(len=64)::filnam,fl5

do i=1,nid
  if (yv(i,1).ne.-99999) then
    tyv(i,1)=yv(i,1)
  else
    print*,"check"
  end if
end do

nid=0
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
do
  nid=nid+1
  read(35,*,iostat=io)cdum4
  if (io.ne.0) exit
  !print*,trim(pedor(i,1)),trim(pedor(i,2))
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
!print*,nid,nm

ALLOCATE(mpos(nm),cno(nm),mafreq1(nm,2),mafreq2(nm)) 
allocate(pedor(nid,4))
allocate(sqgenom(int(nid/4)+1,nm),sex(nid))
allocate(markt(nm),htz(nm))

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


  PRINT*,'*********************************************************************'
  PRINT*,'allele frequency'
  !finding no. marker alleles
  !$OMP PARALLEL DO PRIVATE (yi,xi,w1) 
  do im=1,nm
    x1=0;x2=0;x3=0
    do i=1,nid
      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom(yi,im),xi*2,1),ibits(sqgenom(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        x1=x1+w1
        x2=x2+1
        x3=x3+w1**2
      end if
      !kvm(i,im)=w1
    end do
    mafreq2(im)=x1/x2     !mean*2
    htz(im)=(x3-(x1**2/x2))/(x2) 
  end do
  !$OMP END PARALLEL DO
  print*,'done ****************************************************************' 

  allocate(kvm(nid,nm),xx(nm,nm),xy(nm,1),sol(nm,1),yhat(nid,1))
  !$OMP PARALLEL DO PRIVATE (yi,xi,w1,x3) 
  do im=1,nm
    do i=1,nid
      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom(yi,im),xi*2,1),ibits(sqgenom(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        kvm(i,im)=(w1-mafreq2(im))*(htz(im)**(0.5*-1))
      else
        kvm(i,im)=0
      end if
    end do
  end do
  !$OMP END PARALLEL DO
  !kvm(:,1)=1

  print*,"(X'X)"
  !call dgemm ('N','T',nid,nid,nm,1.0D0,kvm,nid,kvm,nid,0.0D0,xx,nid)
  call dgemm ('T','N',nm,nm,nid,1.0D0,kvm,nid,kvm,nid,0.0D0,xx,nm)
  print*,"(X'X)-1"
  call cholesky_inv (xx,nm,LA)
  print*,"(X'y)"
  call dgemm ('T','N',nm,1,nid,1.0D0,kvm,nid,tyv,nid,0.0D0,xy,nm)
  do i=1,nm
    print*,xy(i,1),tyv(i,1)
  end do
  print*,"(X'X)-1 (X'y)"
  call dgemm ('N','N',nm,1,nm,1.0D0,xx,nm,xy,nm,0.0D0,sol,nm)
  print*,'y-hat'
  call dgemm ('N','N',nid,1,nm,1.0D0,kvm,nid,sol,nm,0.0D0,yhat,nid)

  open (unit=47,file=trim(fl5)//".sol",status='unknown')
  write(47,'(a3,2a14)')'#','BETA','SE'
  do xi=1,nm
    write(47,'(i8,2f14.5)')xi,sol(xi,1),sqrt(xx(xi,xi))
  end do
  close(47)

  !open (unit=48,file=trim(fl5)//".wmat",status='unknown')
  !do xi=1,nid
  !  write(48,'(10000f14.5)')kvm(xi,:)
  !end do
  !close(48)

  !R2 calculations
  sst=0;sse=0
  do i=1,nid
    sst=sst+tyv(i,1)**2
    sse=sse+(tyv(i,1)-yhat(i,1))**2
  end do
  print*,1-sse/sst,1-(sse/(nid-nm-1))/(sst/(nid-1))
  

end subroutine

