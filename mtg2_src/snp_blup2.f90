
subroutine sblup2 (filnam,filnam2,filnam3,flpm,flgpg,fl4,filnam4,sclf)
 
implicit none
integer::nm,nid,sn,dn,i,j,k,zi,zj,im,ma1,wi
integer::ii,idum,nm2,nm3,xi,yi,xi2,yi2
!integer*1,allocatable::igenom(:,:)
real::sclf  !scale factor

character(len=100),allocatable::pedor(:,:)
character(len=50),allocatable::markt(:) !,markt2(:)
real,allocatable::mafreq2(:),mpos(:),htz(:)
!real,allocatable::kvm(:,:),grm(:,:)
double precision,allocatable::kvm(:,:),kvm2(:,:),grm(:,:),snp_blup(:,:)
double precision,allocatable::pm(:,:),vcva(:,:,:),vcve(:,:),snp_vm(:,:)
double precision,allocatable::sblup_r(:,:),kvm3(:,:)
double precision::dx1,sblup_v,sblup_pev

!character(len=1),allocatable::genom(:,:),mafreq1(:,:) !,ran(:)
character(len=1),allocatable::mafreq1(:,:) !,ran(:)
character(len=3)::cdum4,cdum
character(len=3),allocatable::cno(:)

integer::ip,jp,icno,optn

real::rdum(1,1),x1,x2,x3,rij2,wij,wij2,fra,frb,frc,frd,sac,sab,sad,sbc,sbd,scd
real::rij,thtz,ssm,ssm2,v1

character(len=64),allocatable::cdum_t(:) !for each trait
double precision,allocatable::yv(:,:) !,g_blup(:,:)

!sequential storing
integer*1,allocatable::sqgenom(:,:),sex(:) !,igenom(:,:)
integer*1 b1,igen(0:1,0:1),w1,w2
integer::iv1

!progressing report **************************
character*3 bsp

!reading command line *****************
integer::narg,io
character(len=64)::filnam,filnam2,filnam3,flpm,flgpg,filnam4,fl4
logical::file_exist

!CDF
real ( kind = 8 ) bound,mean,p,q,sd,x
integer ( kind = 4 ) status,which


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

allocate(mpos(nm),cno(nm),mafreq1(nm,2),mafreq2(nm)) 
allocate(pedor(nid,4))
allocate(sqgenom(int(nid/4)+1,nm),sex(nid))
allocate(snp_blup(nm,1),htz(nm))
allocate(kvm(nid,nm),kvm2(nm,nid),grm(nid,nid))
allocate(markt(nm))
allocate(cdum_t(1),yv(nid,1))

!if (filnam3=="a") then

open (unit=37,file=trim(filnam3),status='old')      !VgPy file
do i=1,nid
  read(37,*)cdum_t(:)

  !do j=1,tn
  do j=1,1
    if (cdum_t(j).ne.'NA') then
      read(cdum_t(j),*)yv(i,j)
    else
      !yv(i,j)=-99999   !missing
      yv(i,j)=0   ! 0 has no contribution because W' %*% VgPy
    end if
  end do

  !print*,y(i,:)
  !pause
end do
close(37)

!end if


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
  read(35,*)pedor(i,1:4),sex(i),cdum4
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
  if (nm<100 .or. mod(im,int(real(nm)/100))==0) then
    write (*,'(a3,i3)',advance='no')bsp,int(real(im)/real(nm)*100)
    !print*,im
!20  format (A3,I8)
  end if

  do i=1,int((nid-1)/4)+1    !no. ID / 4
    read(38)sqgenom(i,im)
  end do
  !print*,igenom(i,1:20)
  !pause

end do
write (*,'(a22)')'% read plink files'    !closing progress report

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


!reading Ginv **********************
open (unit=71,file=trim(filnam2),status='old')
print*,'Ginv file:',trim(filnam2)
do
  read (71,*,iostat=io)i,j,x1   !order same as in .fam (check in RTMX)
  if (io.ne.0) exit
   grm(i,j)=x1
   grm(j,i)=x1
 end do
close(71)

!indexed genotype coefficients into igen********************************
igen(0,0)=2
igen(0,1)=1
igen(1,1)=0
igen(1,0)=3           !missing *****************************************

PRINT*,'*********************************************************************'
PRINT*,'allele frequency'
open (unit=43,file=trim(filnam)//".freq",status='old') 
do im=1,nm
  read(43,*)cdum4,cdum4,mafreq2(im),htz(im)
end do
close(43)

!finding no. marker alleles
  !g_blup=0
  !$OMP PARALLEL DO PRIVATE (yi,xi,w1)
  do im=1,nm
    do i=1,nid
      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom(yi,im),xi*2,1),ibits(sqgenom(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        kvm(i,im)=(w1-mafreq2(im))*(htz(im)**(0.5*sclf))
      else
        kvm(i,im)=0
      end if
    end do
  end do
  !$OMP END PARALLEL DO
  !print*,kvm
  !print*,grm
  call dgemm ('T','N',nm,nid,nid,1.0D0,kvm,nid,grm,nid,0.0D0,kvm2,nm)
  call dgemm ('N','N',nm,1,nid,1.0D0,kvm2,nm,yv,nid,0.0D0,snp_blup,nm)

PRINT*,'*********************************************************************'
print*,'snp blup done'
  
!if (filnam3=="a") then
if (flpm=="null") then
  open (UNIT=41,FILE=trim(filnam4),STATUS='unknown')
  do i=1,nm
    !write(41,*)markt(i),mafreq1(i,1),snp_blup(i,:)/nm
    write(41,'(a12,a3,f22.16)')markt(i),mafreq1(i,1),snp_blup(i,1)/nm
  end do
  close(41)
else

  allocate(pm(nid,nid),snp_vm(nm,nm)) !check pm dimension should be improved
  allocate(vcva(1,1,1),vcve(1,1),sblup_r(nm,1))  !check dimension
  !reading P mat **********************
  open (unit=72,file=trim(flpm),status='old')
  print*,'P mat file:',trim(flpm)
  do
    read (72,*,iostat=io)i,j,x1   !order same as in .fam (check in RTMX)
    if (io.ne.0) exit
     pm(i,j)=x1
     pm(j,i)=x1
   end do
  close(72)

  inquire(file=fl4,exist=file_exist)
  if (file_exist) then
    open (UNIT=44,FILE=fl4,STATUS='old')
    do i=1,1  !tn   !check # traits
      read(44,*,iostat=io)cdum,vcve(i,i)
      if (io.ne.0) exit
    end do
    do wi=1,1 !mn    !check # random effets
      do xi=1,1 !tn  !check
        read(44,*,iostat=io)cdum,vcva(wi,xi,xi)
        if (io.ne.0) exit
      end do
      do xi=1,1 !tn  !check
        do yi=1,xi-1
          read(44,*,iostat=io)cdum,vcva(wi,xi,yi)
          if (io.ne.0) exit
        end do
      end do
    end do
    close(44)
  else
    print*,'starting value have to be provided for this function >>>'
    pause
  end if

  dx1=vcva(1,1,1)/nm
  call dgemm ('T','N',nm,nid,nid,dx1,kvm,nid,pm,nid,0.0D0,kvm2,nm)
  call dgemm ('N','N',nm,nm,nid,dx1,kvm2,nm,kvm,nid,0.0D0,snp_vm,nm)

  open (UNIT=41,FILE=trim(filnam4),STATUS='unknown')
  !write(41,'(a10,4a13,2a13)')'  SNP_#   ','SNP_EBVs','r2','SE','PEV','Chi','P'
  write(41,'(a10,a7,4a13,2a13)')'SNP','Ref_A','SNP_EBVs','r2','SE','PEV','Chi','P'
    do i=1,nm
      sblup_r(i,1)=snp_vm(i,i)/dx1
      sblup_v=sblup_r(i,1)*dx1
      sblup_pev=(1-sblup_r(i,1))*dx1 
      !CDF functions *****************************************************
      p = huge ( p )
      q = huge ( q )
      sd = 1.0D+00 !degree of greedom for cdfchi
      which=1  !Calculate P and Q from X, MEAN and SD;
      x = (snp_blup(i,1)/nm)**2/sblup_v
      call cdfchi ( which, p, q, x, sd, status, bound )
      !write(41,'(i8,4f14.5,f12.3,e14.5)')i,real(snp_blup(i,1)/nm),real(sblup_r(i,1)),real(sqrt(sblup_v)),real(sblup_pev),real(x),real(q)
      write(41,'(a12,a3,4f14.5,f12.3,e14.5)')markt(i),mafreq1(i,1),real(snp_blup(i,1)/nm),real(sblup_r(i,1)),real(sqrt(sblup_v)),real(sblup_pev),real(x),real(q)
    end do
  close(41)

  deallocate(pm,vcva,vcve,snp_vm,sblup_r)
end if

!ssGWAS ********************************************************************
if (flgpg.ne."null") then    !ssGWAS
  if (flpm.ne."null") then
    print*,"ssGWAS cannot be done with usual SNP BLUP >>> -pm to be removed"
    pause
  endif

  allocate(pm(nid,nid),snp_vm(nm,nm),kvm3(nm,nid)) !check gpg dimension
  allocate(vcva(1,1,1),vcve(1,1),sblup_r(nm,1))  !check dimension
  !reading GPG mat **********************
  open (unit=73,file=trim(flgpg),status='old')
  print*,'GPG mat file:',trim(flgpg)
  do
    read (73,*,iostat=io)i,j,x1   !order same as in .fam (check in RTMX)
    if (io.ne.0) exit
     pm(i,j)=x1
     pm(j,i)=x1
   end do
  close(73)

  inquire(file=fl4,exist=file_exist)
  if (file_exist) then
    open (UNIT=44,FILE=fl4,STATUS='old')
    do i=1,1  !tn   !check # traits
      read(44,*,iostat=io)cdum,vcve(i,i)
      if (io.ne.0) exit
    end do
    do wi=1,1 !mn    !check # random effets
      do xi=1,1 !tn  !check
        read(44,*,iostat=io)cdum,vcva(wi,xi,xi)
        if (io.ne.0) exit
      end do
      do xi=1,1 !tn  !check
        do yi=1,xi-1
          read(44,*,iostat=io)cdum,vcva(wi,xi,yi)
          if (io.ne.0) exit
        end do
      end do
    end do
    close(44)
  else
    print*,'starting value have to be provided for this function >>>'
    pause
  end if

  dx1=vcva(1,1,1)/nm
  call dgemm ('N','N',nm,nid,nid,1.0D0,kvm2,nm,pm,nid,0.0D0,kvm3,nm)
  call dgemm ('N','T',nm,nm,nid,1.0D0,kvm3,nm,kvm2,nm,0.0D0,snp_vm,nm)
  !print*,snp_vm(1,1:5)
  !print*,snp_vm(2,1:5)

  open (UNIT=41,FILE=trim(filnam4),STATUS='unknown')
  write(41,'(a10,a7,4a13,2a13)')'SNP','Ref_A','SNP_EBVs','r2','SE','PEV','Chi','P'
    do i=1,nm
      sblup_r(i,1)=snp_vm(i,i)/dx1/nm**2
      sblup_v=sblup_r(i,1)*dx1
      sblup_pev=(1-sblup_r(i,1))*dx1 
      !CDF functions *****************************************************
      p = huge ( p )
      q = huge ( q )
      sd = 1.0D+00 !degree of greedom for cdfchi
      which=1  !Calculate P and Q from X, MEAN and SD;
      x = (snp_blup(i,1)/nm)**2/sblup_v
      call cdfchi ( which, p, q, x, sd, status, bound )
      write(41,'(a12,a3,4f14.5,f12.3,e14.5)')markt(i),mafreq1(i,1),real(snp_blup(i,1)/nm),real(sblup_r(i,1)),real(sqrt(sblup_v)),real(sblup_pev),real(x),real(q)
    end do
  close(41)

  deallocate(pm,vcva,vcve,snp_vm,sblup_r,kvm3)
endif


end subroutine

