
subroutine sblup (filnam,filnam2,filnam3,flpm,fl4,filnam4,sclf)
 
implicit none
integer::nm,nid,sn,dn,i,j,k,k0,k1,k2,k3,k4,zi,zj,im,ma1
integer::ii,idum,nm2,nm3,xi,yi,xi2,yi2,snpvn
integer,allocatable::flip(:)
!integer*1,allocatable::igenom(:,:)
real::sclf  !scale factor

character(len=100),allocatable::pedor(:,:)
character(len=50),allocatable::markt(:),markt2(:)
real,allocatable::mafreq2(:),mpos(:),htz(:),mhtz(:,:),exv(:,:)
!real,allocatable::kvm(:,:)
double precision,allocatable::kvm(:,:),kvm2(:,:),yv(:,:)
double precision,allocatable::pm(:,:),vcva(:,:,:),vcve(:,:),snp_vm(:,:)
double precision,allocatable::sblup_r(:,:),snp_blup(:,:)
double precision::dx1,sblup_v,sblup_pev

character(len=1),allocatable::genom(:,:),mafreq1(:,:),ran(:)
character(len=3)::cdum4
character(len=3),allocatable::cno(:)

integer::ip,jp,icno,optn,wi

real::rdum(1,1),x1,x2,x3,rij2,wij,wij2,fra,frb,frc,frd,sac,sab,sad,sbc,sbd,scd
real::rij,thtz,ssm,ssm2,v1

character(len=64),allocatable::cdum_t(:) !for each trait
real,allocatable::snpv(:,:),g_blup(:,:)

!sequential storing
integer*1,allocatable::sqgenom(:,:),sex(:) !,igenom(:,:)
integer*1 b1,igen(0:1,0:1),w1,w2
integer::iv1

!progressing report **************************
character*3 bsp

!reading command line *****************
integer::narg,io
character(len=64)::filnam,filnam2,filnam3,filnam4,flpm,fl4
!character(len=128)::filnam,filnam2,filnam3,filnam4

character(len=10000)::string
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

snpvn=1  !20
if (filnam3=="b") then
  open (unit=37,file=trim(filnam2),status='old')
  read(37,'(a)')string
  close(37)

  do i=1,1000
    read(string,*,iostat=io)(cdum4,j=1,i)
    if (io.ne.0) exit
  end do
  snpvn=i-1
  snpvn=snpvn-2
  if (snpvn == 2) then
    print*,'SE of SNP BLUP detected, otherwise check >>>'
  end if
  if (snpvn > 2) then
    print*,'SNP BLUP and SE of SNP BLUP, what else? check >>>'
  end if
  !print*,'# SNP effects:',snpvn
  print*,'NOTE: The same SNPs (order) to be used for bim file and ',trim(filnam2)
  print*,''
end if


ALLOCATE(mpos(nm),cno(nm),mafreq1(nm,2),mafreq2(nm),ran(nm)) 
allocate(pedor(nid,4),flip(nm))
allocate(sqgenom(int(nid/4)+1,nm),sex(nid))
allocate(snp_blup(nm,1),snpv(nm,snpvn),htz(nm))
allocate(g_blup(nid,snpvn))
!allocate(g_blup(nid,1))
allocate(kvm(nid,nm))
!allocate(kvm(nid,1))
allocate(markt(nm),markt2(nm))
allocate(cdum_t(1),yv(nid,1))

if (filnam3=="a") then

open (unit=37,file=trim(filnam2),status='old')      !VgPy file
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

elseif (filnam3=="b") then

open (unit=37,file=trim(filnam2),status='old')      !VgPy file
do i=1,nm
  read(37,*)markt2(i),ran(i),snpv(i,:)
end do
close(37)


end if


open (unit=38,file=trim(filnam)//'.bed',status='old',access='stream',form='unformatted')
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file

flip=0;k0=0;k1=0;k2=0;k3=0;k4=0
do i=1,nm
  read(36,*)cno(i),markt(i),cdum4,mpos(i),mafreq1(i,1),mafreq1(i,2)
  if (filnam3=="b") then
    if (markt(i).ne.markt2(i)) then
      print*,"marker not matched",i,markt(i),markt2(i)
      pause
    end if
    if (mafreq1(i,1)=="G" .and. mafreq1(i,2)=="C") then
      flip(i)=4
      k4=k4+1
      !print*,i,mafreq1(i,1),mafreq1(i,2)
    elseif (mafreq1(i,1)=="C" .and. mafreq1(i,2)=="G") then
      flip(i)=4
      k4=k4+1
      !print*,i,mafreq1(i,1),mafreq1(i,2)
    elseif (mafreq1(i,1)=="T" .and. mafreq1(i,2)=="A") then
      flip(i)=4
      k4=k4+1
      !print*,i,mafreq1(i,1),mafreq1(i,2)
    elseif (mafreq1(i,1)=="A" .and. mafreq1(i,2)=="T") then
      flip(i)=4
      k4=k4+1
      !print*,i,mafreq1(i,1),mafreq1(i,2)
    else
      if (mafreq1(i,1)==ran(i)) then
        k0=k0+1
      elseif (mafreq1(i,1).ne.ran(i)) then
        !print*,"reference allele not matched"
        !pause
        if (mafreq1(i,2)==ran(i)) then
          flip(i)=1
          k1=k1+1
        else
          if (ran(i)=="G") then
            if (mafreq1(i,1)=="C") then
              flip(i)=2
              k2=k2+1
            elseif (mafreq1(i,2)=="C") then
              flip(i)=3
              k3=k3+1
            end if
          elseif (ran(i)=="C") then
            if (mafreq1(i,1)=="G") then
              flip(i)=2
              k2=k2+1
            elseif (mafreq1(i,2)=="G") then
              flip(i)=3
              k3=k3+1
            end if
          elseif (ran(i)=="T") then
            if (mafreq1(i,1)=="A") then
              flip(i)=2
              k2=k2+1
            elseif (mafreq1(i,2)=="A") then
              flip(i)=3
              k3=k3+1
            end if
          elseif (ran(i)=="A") then
            if (mafreq1(i,1)=="T") then
              flip(i)=2
              k2=k2+1
            elseif (mafreq1(i,2)=="T") then
              flip(i)=3
              k3=k3+1
            end if
          end if
        end if
      end if
    end if
  end if
  !print*,mafreq1(i,:)
end do

if (filnam3=="b") then
  print*,"1) reference matched              :",k0
  print*,"2) reference flip SNP # (flipped) :",k1
  print*,"3) strand flip SNP # (flipped)    :",k2
  print*,"4) strand flip + reference flip   :",k3
  print*,"5) ambiguous SNPs                 :",k4
  print*,"NOTE: Carefully check strand/reference alleles if 2)-5) is non zero"
  print*,''
end if

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
open (unit=43,file=trim(filnam)//".freq",status='old') 

!finding no. marker alleles
  g_blup=0
  do im=1,nm

    if (nm<=100 .or. mod(im,int(nm/100))==0) then
      write (*,'(a3,i3)',advance='no')bsp,int((.1*im)/(.1*nm)*100)
    end if

    read(43,*)cdum4,cdum4,mafreq2(im),htz(im)
    !htz(im)=2*(mafreq2(im)/2)*(1-mafreq2(im)/2)

    x1=0;x2=0
    do i=1,nid

      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom(yi,im),xi*2,1),ibits(sqgenom(yi,im),xi*2+1,1))
      if (w1.ne.3) then
      !if (w1.ne.3.and.yv(i,1).ne.0) then
        x1=x1+w1
        x2=x2+1
      end if
      kvm(i,im)=w1
      !kvm(i,1)=w1

      !reference flip for optin b *******************
      if (filnam3=="b") then
        if (flip(im)==1 .or. flip(im)==3) then
          !kvm(i,1)=2-kvm(i,1)
          kvm(i,im)=2-kvm(i,im)
          !print*,mafreq2(im)
          mafreq2(im)=(1-mafreq2(im)/2)*2
          !print*,mafreq2(im)
          !pause
        end if
      end if !***************************************
    end do

    do i=1,nid
      !if (kvm(i,1)==3) then
      if (kvm(i,im)==3) then
        !kvm(i,1)=0
        kvm(i,im)=0
      else
        !kvm(i,1)=(kvm(i,1)-mafreq2(im))/htz(im)**0.5
        !kvm(i,im)=(kvm(i,im)-mafreq2(im))/htz(im)**0.5
        kvm(i,im)=(kvm(i,im)-mafreq2(im))*(htz(im)**(0.5*sclf))
      end if

      if (filnam3=="b") then
        !do ii=1,snpvn
        do ii=1,1
          !g_blup(i,ii)=g_blup(i,ii)+kvm(i,1)*snpv(im,ii)
          g_blup(i,ii)=g_blup(i,ii)+kvm(i,im)*snpv(im,ii)
        end do
      end if

    end do

    !if (filnam3=="a") then
    !  rdum=matmul(transpose(kvm),yv)
    !  snp_blup(im,1)=rdum(1,1)
    !end if
  end do
  write (*,'(a6)')'% done'    !closing progress report
  !print*,mafreq2(1) !,thtz,ssm

close(43)

PRINT*,'*********************************************************************'
print*,'snp blup '

  
if (filnam3=="a") then

  !print*,kvm(1,1:5)
  !print*,kvm(2,1:5)
  !print*,kvm(3,1:5)
  call dgemm ('T','N',nm,1,nid,1.0D0,kvm,nid,yv,nid,0.0D0,snp_blup,nm)
  !print*,snp_blup(1:5,1)

  if (flpm=='null') then

  !open (UNIT=41,FILE=trim(filnam)//'.snpv',STATUS='unknown')
  open (UNIT=41,FILE=trim(filnam4),STATUS='unknown')
  do i=1,nm
    !write(41,*)markt(i),mafreq1(i,1),snp_blup(i,:)/nm
    write(41,'(a12,a3,f22.16)')markt(i),mafreq1(i,1),snp_blup(i,1)/nm
  end do
  close(41)

  else

  allocate(pm(nid,nid),snp_vm(nm,nm),kvm2(nm,nid)) !check pm dimension should be improved
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
      read(44,*,iostat=io)cdum4,vcve(i,i)
      if (io.ne.0) exit
    end do
    do wi=1,1 !mn    !check # random effets
      do xi=1,1 !tn  !check
        read(44,*,iostat=io)cdum4,vcva(wi,xi,xi)
        if (io.ne.0) exit
      end do
      do xi=1,1 !tn  !check
        do yi=1,xi-1
          read(44,*,iostat=io)cdum4,vcva(wi,xi,yi)
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
      !write(41,'(a12,a3,4f14.5,f12.3,e14.5)')markt(i),mafreq1(i,1),real(snp_blup(i,1)/nm),real(sblup_r(i,1)),real(sqrt(sblup_v)),real(sblup_pev),real(x),real(q)
      write(41,'(a40,a3,4f14.5,f12.3,e14.5)')markt(i),mafreq1(i,1),real(snp_blup(i,1)/nm),real(sblup_r(i,1)),real(sqrt(sblup_v)),real(sblup_pev),real(x),real(q)
    end do
  close(41)

  deallocate(pm,vcva,vcve,snp_vm,sblup_r)
  end if


elseif (filnam3=="b") then
  if (snpvn==1) then
    open (UNIT=41,FILE=trim(filnam4),STATUS='unknown')
    do i=1,nid
      !write(41,'(i8,f22.16)')i,real(g_blup(i,:))
      write(41,'(i8,f22.16)')i,real(g_blup(i,1))
    end do
    close(41)
  elseif (snpvn==2) then    !reliability 
    do i=1,nid
        x1=0
      do im=1,nm
        x1=x1+kvm(i,im)**2*snpv(im,2)**2  !WDW', D is diagonal with var(SNP_BLUP)
      end do
      g_blup(i,2)=x1
    end do     
    open (UNIT=41,FILE=trim(filnam4),STATUS='unknown')
    do i=1,nid
      write(41,'(i8,2f22.16)')i,real(g_blup(i,1)),real(g_blup(i,2))**0.5
    end do
    close(41)
  end if

end if


end subroutine

