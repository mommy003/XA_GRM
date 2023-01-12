subroutine hrtmx_2 (hrtmx,hrtmx2,grm,pedn,n,mn,fam,tune,blend,method,f41,filnam,f42)

!*********************************************************************
!* To get H matrix for SSGBLUP
!* Authors: S. Hong Lee, University of South Australia (C) 2017-2021
!* (C) 2021 NIAS, RDA - All Rights Reserved under GNU LGPL  
!* Collaborative project (2021-2022) 
!*********************************************************************

implicit none

integer::n,pedn,i,j,w,k,kk,k1,k2,fn,mn,zi,mo,md,tune,im,yi,xi,method
character(len=64)::hrtmx,hrtmx2,f41,filnam,f42
character(len=100)::ped(pedn,4),idx_mem(pedn*5),fam(n,4)  !check
real::grm(mn,n,n),nrm(pedn,pedn),r1,r2,sclf
integer::idx_ped(pedn*5,4),idx_pedn,idx_fam2(n),idx_fam22(n) !x5 check
integer::idx(pedn),idx_fam(n),idx_01(pedn),idx_total(pedn),idx_total2(pedn)
real,allocatable::idx_nrm(:,:)

double precision::A11((pedn-n),(pedn-n)),A21(n,(pedn-n)),GA22(n,n)
double precision::A22(n,n)
integer,allocatable::idx_ped2(:,:)
double precision::x1,x2,x3

double precision::A22A21(n,(pedn-n)),A12A22((pedn-n),n)
double precision::GA22A21(n,(pedn-n)),G22(n,n)
double precision::H11_tmp((pedn-n),n),H11_tmp2((pedn-n),(pedn-n))
double precision::H11_tmp4((pedn-n),n),H11_tmp3((pedn-n),(pedn-n))

double precision,allocatable::wm(:,:),wm_imp(:,:)
real,allocatable::wm_tot(:,:),mafreq2(:,:),htz(:,:)

real::mdiagg,mdiaga,moffg,moffa,av,bv,blend

!plink file for Z1 imputation
integer::nid,nm,io,eoro
character(len=3)::cdum4

!sequential storing
integer*1,allocatable::sqgenom(:,:),sex(:)
integer*1 b1,igen(0:1,0:1),w1,w2


!print*,pedn,n,mn
open (unit=48,file=hrtmx,status='old') 
do i=1,pedn
  read(48,*)ped(i,:)
end do
close(48)

idx_pedn=pedn
idx_mem(1:pedn)=ped(:,2)

!indexing
!!$OMP PARALLEL DO PRIVATE(zi,i)
do zi=1,pedn
  do i=1,idx_pedn
    if (idx_mem(i)==ped(zi,2)) then
      idx(zi)=i
      goto 50
    end if
  end do
50 continue
end do
!!$OMP END PARALLEL DO
print*,'indexing completed ****************************************'


allocate(idx_nrm(idx_pedn,idx_pedn))  !,idx_ped2(idx_pedn,4))
!idx_ped2=idx_ped(1:idx_pedn,:)

!call Amat(idx_nrm,idx_ped2,idx_pedn)
open (unit=49,file=hrtmx2,status='old') 
print*,'NRM file:',trim(hrtmx2)
do
  read(49,*,iostat=io)i,j,x1
  if (io.ne.0) exit
  idx_nrm(i,j)=x1
  idx_nrm(j,i)=x1
end do
close(49)

print*,'NRM completed *********************************************'

r1=0;r2=0;mo=0
do zi=1,pedn
  !do j=1,pedn
  do j=1,zi
    nrm(zi,j)=idx_nrm(idx(zi),idx(j))
    nrm(j,zi)=nrm(zi,j)
    if (zi.ne.j) then
      mo=mo+1
      r2=r2+nrm(zi,j)
    end if
    !print*,nrm(zi,j)
  end do
  r1=r1+nrm(zi,zi)
end do
!print*,r1/pedn,r2/mo,r2/(pedn**2/2-pedn/2)

deallocate(idx_nrm) !,idx_ped2)

if (mn > 0) then  !matching genotype info
  if (mn > 1) then
    print*,'# GRM should not be more than 1'
    pause
  end if

  do zi=1,n
    !print*,zi
    do i=1,pedn
      !print*,i
      if (fam(zi,2) == ped(i,2)) goto 60
    end do
    print*,'no matched ID in pedigree file >>> check',fam(zi,2)
    pause
60  continue
    idx_fam(zi)=i    !ped # for the fam order
  end do
  
  idx_01=0
  do i=1,n
    !print*,grm(1,i,i)
    !print*,idx_fam(i)
    idx_01(idx_fam(i))=1
  end do
  !print*,idx_01

  k=0
  do i=1,pedn
    if (idx_01(i)==0) then
      k=k+1
      idx_total(k)=i
      idx_total2(i)=k        !for getting back to the original order
    end if
  end do
  k1=k
  do i=1,pedn
    if (idx_01(i)==1) then
      k=k+1
      idx_total(k)=i
      idx_total2(i)=k        !for getting back to the original order
    end if
  end do
  k2=k-k1
  !print*,idx_total
  !pause
  !print*,idx_total2

  if (k.ne.pedn .or. k1.ne.pedn-n .or. k2.ne.n) then
    print*,'number not mached',k,k1,k2,pedn,n
    pause
  end if
  print*,'total pedigree #    :',k
  print*,'genotyped samples # :',k2

  !GRM order index for G22
  do i=1,n
    do j=1,n
      if (ped(idx_total(k1+i),2)==fam(j,2)) then
        idx_fam2(i)=j    !fam file # for genotyped ID in ped order
        idx_fam22(j)=i    !for getting back to original fam order
        goto 70
      end if
    end do
    print*,'something wrong >>> check'
70 continue
  end do
  !pause
  !print*,idx_fam2
  print*,'GRM order check ********************************************'

  !A11
  do i=1,k1
    do j=1,k1
      A11(i,j)=nrm(idx_total(i),idx_total(j))
    end do
  end do
  !A21
  do i=1,k2
    do j=1,k1
      A21(i,j)=nrm(idx_total(k1+i),idx_total(j))
    end do
  end do
  !A22
  do i=1,k2
    do j=1,k2
      A22(i,j)=nrm(idx_total(k1+i),idx_total(k1+j))
    end do
  end do
  
  !G22=grm(1,:,:)
  do i=1,k2
    do j=1,k2
      G22(i,j)=grm(1,idx_fam2(i),idx_fam2(j))
    end do
  end do
  print*,'A11, A21, A22 and G22 parts assigned ***********************'

  if (blend.ne.1) then
    print*,'blend G = G*',blend,'+ A22*',1-blend
    do i=1,k2
      do j=1,i
        G22(i,j)=G22(i,j)*blend+A22(i,j)*(1-blend)
        G22(j,i)=G22(i,j)
      end do
    end do
  end if  

  if (tune==1 .or. tune==2) then
    md=0;mo=0
    mdiagg=0;mdiaga=0;moffg=0;moffa=0
    do i=1,k2
      do j=1,i
        if (i.ne.j) then
          mo=mo+1
          moffg=moffg+G22(i,j)
          moffa=moffa+A22(i,j)
        end if
      end do
      mdiagg=mdiagg+G22(i,i)
      mdiaga=mdiaga+A22(i,i)
    end do
    mdiagg=mdiagg/k2
    mdiaga=mdiaga/k2
    moffg=moffg/mo
    moffa=moffa/mo

    bv=(mdiaga-moffa)/(mdiagg-moffg)
    av=((k2**2-k2)*moffa+k2*mdiaga)/k2**2 - ((k2**2-k2)*moffg+k2*mdiagg)/k2**2

    print*,'dG-oG',mdiagg-moffg,'dA-oA',mdiaga-moffa
    print*,'G all',((k2**2-k2)*moffg+k2*mdiagg)/k2**2,'A all',((k2**2-k2)*moffa+k2*mdiaga)/k2**2

    if (tune==1) then
      print*,'tuning G - following Method 2 in blupf90' 
      print*,'final G =',av,'+ (blended) G*',bv
      do i=1,k2
        do j=1,i
          G22(i,j)=av+G22(i,j)*bv
          G22(j,i)=G22(i,j)
        end do
      end do
    elseif (tune==2) then
      print*,'tuning G - following Method 4 in blupf90' 
      print*,'final G =',av,'+ (blended) G*1 '
      do i=1,k2
        do j=1,i
          G22(i,j)=av+G22(i,j)
          G22(j,i)=G22(i,j)
        end do
      end do

    end if
  end if

  !if (tune==1.or.blend.ne.1) then
  if (f42.ne.'null') then
    !open (unit=51,file=trim(f41)//'.grm',status='unknown')
    open (unit=51,file=trim(f42),status='unknown')
    do zi=1,n
      do j=1,zi
        !if (zi==j) then
        !  write(51,'(i8,i8,f14.8)')zi,j,grm(1,zi,j)+(mdiaga-mdiagg)
        !else
        !  write(51,'(i8,i8,f14.8)')zi,j,grm(1,zi,j)+(moffa-moffg)
        !end if
        write(51,'(i8,i8,f14.8)')zi,j,G22(idx_fam22(zi),idx_fam22(j))
      end do
    end do
    close(51)
    print*,'blended or/and tuned grm stored ******************************'
  end if

  if (method==1) then
    !G-A22
    GA22=G22-A22
    !print*,GA22

    !A22-1
    call cholesky_inv (A22,k2,x1)
    !print*,A22
    print*,'A22 inverted ***********************************************'

    call dgemm ('N','N',n,(pedn-n),n,1.0D0,A22,n,A21,n,0.0D0,A22A21,n)
    do i=1,(pedn-n)
      do j=1,n
        A12A22(i,j)=A22A21(j,i)
      end do
    end do
    print*,'A12A22 done ************************************************'
    !print*,A12A22
    !pause
    call dgemm ('N','N',(pedn-n),n,n,1.0D0,A12A22,(pedn-n),GA22,n,0.0D0,H11_tmp,(pedn-n))
    call dgemm ('N','N',(pedn-n),(pedn-n),n,1.0D0,H11_tmp,(pedn-n),A22A21,n,0.0D0,H11_tmp2,(pedn-n))
    print*,'A11 part done **********************************************'

    call dgemm ('N','N',n,(pedn-n),n,1.0D0,G22,n,A22A21,n,0.0D0,GA22A21,n)
    !print*,H11_tmp 
    print*,'A21 part done **********************************************'

  elseif (method==2) then
    print*,'method 2'
    H11_tmp2=0
    GA22A21=A21
  end if

  nrm(1:k1,1:k1)=A11+H11_tmp2
  !print*,nrm(3,:)
  do i=1,(pedn-n)
    do j=1,n
      nrm(k1+j,i)=GA22A21(j,i)
      nrm(i,k1+j)=GA22A21(j,i)
    end do
  end do
  !print*,nrm(3,:)
  nrm((k1+1):(k1+k2),(k1+1):(k1+k2))=G22

  !open (unit=50,file=trim(hrtmx)//'.hrtmx',status='unknown')
  open (unit=50,file=trim(f41),status='unknown')
  do zi=1,pedn
    do j=1,zi
      if (nrm(idx_total2(zi),idx_total2(j)) .ne. 0) then
        write(50,'(i8,i8,f14.8)')zi,j,nrm(idx_total2(zi),idx_total2(j))
        !print*,idx_total2(zi),idx_total2(j)
      end if
    end do
  end do
  close(50)
  print*,'H matrix stored ********************************************'

  if (filnam .ne. 'null') then  !plink file for Z1 imputation

nid=0
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
do
  nid=nid+1
  read(35,*,iostat=io)cdum4
  if (io.ne.0) exit
end do
close(35)
nid=nid-1
if (n .ne. nid) then
  print*,'no. ID not matched >>> check',nid,n
end if

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

allocate(sqgenom(int(nid/4)+1,nm),wm(nid,nm),wm_imp((pedn-nid),nm))
allocate(mafreq2(nm,1),htz(nm,1),wm_tot(pedn,nm))

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
    read(38)sqgenom(i,im)
  end do
  !print*,igenom(i,1:20)
  !pause
end do
write (*,'(a22)')'100% read plink files'    !closing progress report
close(38)


!indexed genotype coefficients into igen********************************
igen(0,0)=2
igen(0,1)=1
igen(1,1)=0
igen(1,0)=3           !missing *****************************************

PRINT*,''

eoro=2    !check
do zi=1,1
  PRINT*,'allele frequency' ! for pop',zi
  !$OMP PARALLEL DO PRIVATE (yi,xi,w1) 
  do im=1,nm
    x1=0;x2=0;x3=0
    do i=1,nid  !pedn
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
    !print*,mafreq2(im,zi),htz(im,zi)
  end do
  !$OMP END PARALLEL DO
end do

  sclf=-1       !check
  !allocate(wm(nid,nm))
  !$OMP PARALLEL DO PRIVATE (yi,xi,w1,x3) 
  do im=1,nm
    do i=1,nid
      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom(yi,im),xi*2,1),ibits(sqgenom(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        wm(i,im)=(w1-mafreq2(im,1))*(htz(im,1)**(0.5*sclf))
      else
        wm(i,im)=0
      end if
    end do
  end do
  !print*,wm(1,:)
  !print*,wm(2,:)
  !print*,wm(3,:)
  !Z1 imputation 
  call dgemm ('N','N',(pedn-n),nm,n,1.0D0,A12A22,(pedn-n),wm,n,0.0D0,wm_imp,(pedn-n))

  wm_tot(1:k1,:)=wm_imp(1:k1,:)
  wm_tot((k1+1):pedn,:)=wm(1:n,:)
  !print*,wm_imp

  open (unit=52,file=trim(f41)//'.wmat',status='unknown')
  do zi=1,pedn
    write(52,'(100000000f14.8)')wm_tot(idx_total2(zi),:)
  end do
  close(52)
  print*,'Imputed W1 + W2 (standardised) matrix stored in ',trim(f41)//'.wmat'

  end if ! if (filnam .ne. 'null') 


else ! if (mn ==0)

  !open (unit=50,file=trim(hrtmx)//'.hrtmx',status='unknown')
  open (unit=50,file=trim(f41),status='unknown')
  do zi=1,pedn
  !do zi=1,n
    do j=1,zi
      if (nrm(zi,j) .ne. 0) then
        !write(50,*)zi,j,nrm(zi,j)
        write(50,'(i8,i8,f14.8)')zi,j,nrm(zi,j)
      end if
      !if (A22(zi,j) .ne. 0) then
      !  write(50,*)zi,j,A22(zi,j)
      !end if
    end do
  end do
  close(50)
  print*,'A matrix stored ********************************************'

end if !mn > 0



end subroutine
