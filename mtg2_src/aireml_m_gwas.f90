
subroutine aireml_m_gwas (mbin,yidx,yv,rn,rnm,trn,pedn,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,fl24,xmm,nit,conv,wt_res,filnam)        

!***********************************************************************
!aireml_b: multivariate analysis especially for case control 
!S. Hong Lee

!mbin : matrices bin
!obs  : observed ID
!yv   : phenotypes
!rn   : no. phenotypes
!pedn : no. pedigree (ID) should be = rn
!fn   : no. fixed effects
!mn   : no. matrices (random effects)
!***********************************************************************

implicit none

INTEGER::n,nped,ix,iz,tn,nb,nbc,nbc2,xi,yi,ri,zi,vi,ui,wi,wj,nit
integer::rn,pedn,mn,trn,trnx,trny,tfn,trnv,trnu,itit
integer::yidx(tn,rn),rnm(tn),fnm(tn)
integer::rn1,rn2,yidx1(rn),yidx2(rn)
!double precision::mbin(mn,pedn,pedn),yv(rn,2),tyv(rn*2,1)
real::mbin(mn,pedn,pedn),yv(rn,tn),tyv(trn,1),sv_fact
double precision::resw(trn),resw2(pedn,tn)

!double precision::qq(mn,rn,rn),qq2(mn,rn,rn)
!real::qq(mn,rn,rn),qq2(mn,rn,rn)

double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),conv
double precision::vcva(mn,tn,tn),vcve(tn,tn),vvc(tn)
double precision::fix_vcva(mn,tn,tn),fix_vcve(tn,tn)
double precision::vcva_rsv(mn,tn,tn),vcve_rsv(tn,tn)
double precision::blup_ebv(pedn,tn,mn),blup_py(pedn,tn,mn),beta(tfn,1) !,cor
double precision::blup_r(pedn,tn,mn)

integer::i,j,m,k,l,io
double precision::sum_h2,sum_v(tn)


double precision::LC,LA,LG,ypy,x1,x2,x3,x10,x11,v1,v2,v10,v11
double precision::y1,y2,y3,y10,y11,z1,z2,z3,z10,z11
double precision::LKHP,MLKH,mva(mn),mva2(mn),mve,mve2,mcov(mn)
double precision::h2(mn),h2n(mn),tr1,tr2(mn),tr3

double precision::xm(trn,tfn),xmm(tn,pedn,10000)

double precision::xvmx(tfn,tfn),xvm(tfn,trn) 
double precision::xvm2(tfn,trn),xvm3(trn,tfn) 
double precision::pm(trn,trn),py(trn,1),xb(trn,1)

double precision::v_tmp(trn),v2_tmp(trn)

!double precision::aim(2+mn*3,2+mn*3),aiminv(2+mn*3,2+mn*3),up(2+mn*3,1),dldv(2+mn*3,1)
!double precision::sdm(2+mn*3,2+mn*3),sdmI(2+mn*3,2+mn*3)
double precision::aim(tn+mn*(tn+(tn**2-tn)/2),tn+mn*(tn+(tn**2-tn)/2))
double precision::sdmI(tn+mn*(tn+(tn**2-tn)/2),tn+mn*(tn+(tn**2-tn)/2))
double precision::up(tn+mn*(tn+(tn**2-tn)/2),1),dldv(tn+mn*(tn+(tn**2-tn)/2),1)

!for BLUP ebv
!double precision::vmat(mn,pedn,rn)

!time measure **********************************************************
INTEGER::now(8)
CHARACTER*8 date
CHARACTER*20 time,zone
double precision:: timef, t1_real, t2_real, elapsed_real_secs
double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
!***********************************************************************

character(len=7)::cdum
character(len=64)::fl4,fl24,fl5,fl11,fl11r,wt_res,cdum_t(tn),filnam
logical::file_exist

!sequential storing
integer*1,allocatable::sqgenom(:,:),sex(:) !,igenom(:,:)
integer*1 b1,igen(0:1,0:1),w1,w2

integer::nid,nsnp,im,gwasi
integer,allocatable::snpv(:),nmiss(:)   
character(len=20),allocatable::mnam(:),markt(:),cmk(:)
real,allocatable::mpos(:),raf(:)

character(len=1),allocatable::mafreq1(:,:)
character(len=3)::cdum4
character(len=3),allocatable::cno(:)

!CDF
real ( kind = 8 ) bound,mean,p,q,sd,x
integer ( kind = 4 ) status,which


!progressing report **************************
character*3 bsp
!character*1 creturn



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


nsnp=0
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file
do
  nsnp=nsnp+1
  read(36,*,iostat=io)cdum4
  if (io.ne.0) exit
end do
close(36)
nsnp=nsnp-1
print*,'no. marker: ',nsnp

ALLOCATE(mpos(nsnp),cno(nsnp),mnam(nsnp),mafreq1(nsnp,2))
!allocate(pedor(nid,4))
allocate(sqgenom(int(nid/4)+1,nsnp),sex(nid))
allocate(markt(nsnp),snpv(nid),raf(nsnp),nmiss(nsnp))


open (unit=38,file=trim(filnam)//'.bed',status='old',access='stream',form='unformatted')
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file


do i=1,nsnp
  read(36,*)cno(i),markt(i),cdum4,mpos(i),mafreq1(i,1),mafreq1(i,2)
  !print*,mafreq1(i,:)
end do

!progressing report***************************************************
bsp(1:1)=char(8); bsp(2:2)=char(8); bsp(3:3)=char(8)

!!for fam file (ped) ***************************************************
!do i=1,nid
!  read(35,*)pedor(i,:),sex(i),cdum4
!end do !***************************************************************

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

do im=1,nsnp
  if (nsnp<100 .or. mod(im,int(nsnp/100))==0) then
    write (*,'(a3,i3)',advance='no')bsp,int((.1*im)/(.1*nsnp)*100)
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

!!check duplications in ped file **********************************
!do i=1,nid
!  do j=i+1,nid
!    if (pedor(i,1)==pedor(j,1).and.pedor(i,2)==pedor(j,2)) then
!      print*,'ID duplications in .ped: ',trim(pedor(i,1)),trim(pedor(i,2))
!      pause
!    end if
!  end do
!end do !**********************************************************


!indexed genotype coefficients into igen********************************
igen(0,0)=2
igen(0,1)=1
igen(1,1)=0
igen(1,0)=3           !missing *****************************************




xm=0;trn=0;tfn=0
do xi=1,tn
  tfn=tfn+fnm(xi)
  do i=1,rn
    !print*,xm1(i,:)
    if (yv(i,xi).ne.-99999) then
      trn=trn+1
      tyv(trn,1)=yv(i,xi)

      do j=1,fnm(xi)
        xm(trn,tfn-fnm(xi)+j)=xmm(xi,i,j)
      end do
      !print*,xm(rn1,1:fn1),xm1(i,:),rn1,i
      !pause
    end if
  end do
end do
tfn=tfn+tn
print*,tfn

print*,'*** number of records used ***'
do i=1,tn
  print*,'trait',i,':',rnm(i)
end do
print*,''

if (wt_res.ne."null") then
  open (unit=50,file=wt_res,status='old')
  do i=1,pedn
    read(50,*)cdum,cdum,cdum_t(:)
    do j=1,tn
      if (cdum_t(j).ne.'NA') then
        read(cdum_t(j),*)resw2(i,j)
      else
        resw2(i,j)=-99999   !missing
      end if
      !print*,resw2(i,:)
    end do
  end do
  close(50)
  !pause

  k=0
  do j=1,tn
    do i=1,pedn
      if (resw2(i,j).ne.-99999) then
        k=k+1
        resw(k)=1/resw2(i,j)
        !print*,resw(k)
      end if
    end do
    if (k .ne. sum(rnm(1:j))) then
      print*,'# residual weights is not matched with # records'
      pause
    end if
  end do

else      !if null
  resw=1.0D0
end if




open (unit=47,file=trim(fl5)//".gwas",status='unknown')
write(47,'(a3,a17,a3,a8,2a14,a12,a14)')'CHR','SNP','A1','N','BETA','SE','CHI','P'
do gwasi=1,nsnp
  !print'(a4,i8,5x,f5.2,a12)','SNP #',gwasi,real(gwasi)/real(nsnp)*100,"% progressed"
  print'(a4,i8,5x,f5.2,a12)','SNP #',gwasi,real(gwasi)/real(nsnp)*100,"% progressed ***"
  !creturn=achar(13)
  !write (*,'(a,a4,i8,5x,f5.2,a12)',advance='no')creturn,'SNP #',gwasi,real(gwasi)/real(nsnp)*100,"% progressed ***"
  x1=0;x2=0
  do i=1,nid
    yi=int((i-1)/4)+1
    xi=mod((i-1),4)
    w1=igen(ibits(sqgenom(yi,gwasi),xi*2,1),ibits(sqgenom(yi,gwasi),xi*2+1,1))
    snpv(i)=w1
    if (w1.ne.3) then
      x1=x1+w1; x2=x2+1
    end if
  end do
  raf(gwasi)=(x1/x2)*0.5
  nmiss(gwasi)=x2
  do i=1,nid
    if (snpv(i)==3) snpv(i)=x1/x2
  end do

trn=0
do xi=1,tn
  do i=1,rn
    if (yv(i,xi).ne.-99999) then
      trn=trn+1
      xm(trn,tfn-tn+xi)=snpv(yidx(xi,trn))
      !print*,xm(rn1,1:fn1),xm1(i,:),rn1,i
      !pause
    end if
  end do
end do


  ! y = Xb + Zu + e ---- fitting effects without QTL

  sv_fact=0
  !arbitary starting value (half of the total variance)
  vcve=0;vcva=0
  do xi=1,tn
    x1=0;x2=0
    do i=1,rn
      if (yv(i,xi).ne.-99999) then
        x1=x1+yv(i,xi)
        x2=x2+yv(i,xi)**2
      end if
    end do
    vcve(xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))/2
    vcva(:,xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))/(2*mn)
  end do

  
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
  end if

  !print*,'fl24:',fl24
  inquire(file=fl24,exist=file_exist)
  if (file_exist) then
    open (UNIT=51,FILE=fl24,STATUS='old')
    do i=1,tn
      read(51,*,iostat=io)cdum,fix_vcve(i,i)
      if (io.ne.0) exit
    end do
    do wi=1,mn
      do xi=1,tn
        read(51,*,iostat=io)cdum,fix_vcva(wi,xi,xi)
        if (io.ne.0) exit
      end do
      do xi=1,tn
        do yi=1,xi-1
          read(51,*,iostat=io)cdum,fix_vcva(wi,xi,yi)
          if (io.ne.0) exit
        end do
      end do
    end do
    close(51)
  else
    fix_vcve=-99999
    fix_vcva=-99999
  end if


110 continue

  itit=1
  ! Start iteration
  LKH=-1000000000
  MLKH=-1000000000


  ! V matrix    V= ZAZ'Va+Z2GZ2'Vpe+IVe*************************
  !do zi=1,200
  !do zi=1,nit
  do zi=1,10000000
    LKHP=LKH
  !LKH and AIREML (without sparce technique)
  ! V matrix    V = ZAZ'Va + IVe *************************

  call cpu_time (t2_cpu)

  !fixing some parameters
    do i=1,tn
      if (fix_vcve(i,i) .ne. -99999) vcve(i,i)=fix_vcve(i,i)
    end do
    do wi=1,mn
      do xi=1,tn
        if (fix_vcva(wi,xi,xi).ne.-99999) vcva(wi,xi,xi)=fix_vcva(wi,xi,xi)
      end do
      do xi=1,tn
        do yi=1,xi-1
          if (fix_vcva(wi,xi,yi).ne.-99999) vcva(wi,xi,yi)=fix_vcva(wi,xi,yi)
        end do
      end do
    end do
  !**************************
  
  pm=0;trnx=0
  do xi=1,tn
    trnx=trnx+rnm(xi)
    trny=0
    do yi=1,xi
      trny=trny+rnm(yi)
      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do i=1,rnm(xi)
        do j=1,rnm(yi)
          do k=1,mn
            pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)+mbin(k,yidx(xi,i),yidx(yi,j))*vcva(k,xi,yi)
            !print*,pm(trnx-rnm(xi)+i,trny-rnm(yi)+j),trnx-rnm(xi)+i,trny-rnm(yi)+j,mbin(k,yidx(xi,i),yidx(yi,j))
            !pause
          end do
          if (xi.ne.yi) then
            pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)
          end if
        end do
        if (xi==yi) then
          !pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)+vcve(xi,xi)
          pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)+vcve(xi,xi)*resw(trnx-rnm(xi)+i)
        end if
      end do
      !$OMP END PARALLEL DO

    end do
  end do

  call cpu_time (t2_cpu)
  call cholesky_inv (pm,trn,LA)
  call cpu_time (t1_cpu)
  !print*,'V inverse done - time:',real(t1_cpu-t2_cpu) !,LA

    ! P = (VI)-(VI X (X'VI X)I' X' VI)******************************
    !xvm=MATMUL(transpose(xm),pm)

    !$OMP PARALLEL DO PRIVATE (v1)
    do i=1,tfn
      do j=1,trn
        v1=0
        do k=1,trn
          v1=v1+xm(k,i)*pm(k,j)    !xm*pm
        end do
        xvm(i,j)=v1
      end do
    end do
    !$OMP END PARALLEL DO

    !xvmx=MATMUL(xvm,xm)
    !$OMP PARALLEL DO PRIVATE (v1)
    do i=1,tfn
      do j=1,tfn
        v1=0
        do k=1,trn
          v1=v1+xvm(i,k)*xm(k,j)    !xm*pm
        end do
        xvmx(i,j)=v1
      end do
    end do
    !$OMP END PARALLEL DO

    !call cholesky_inv (xvmx,tfn,xvmxI,LG)
    call cholesky_inv (xvmx,tfn,LG)
    !BLUE = (X'VI X)I XVI y **********************
    !beta=matmul(matmul(matmul(xvmx,transpose(xm)),pm),tyv)

    !$OMP PARALLEL DO PRIVATE (v1)
    do i=1,tfn
      do j=1,trn
        v1=0
        do k=1,tfn
          v1=v1+xvmx(i,k)*xm(j,k)    !xvmx*xm
        end do
        xvm2(i,j)=v1
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE (v1,v2)
    do i=1,tfn
      v2=0
      do j=1,trn
        v1=0
        do k=1,trn
          v1=v1+xvm2(i,k)*pm(k,j)    !*pm
        end do
        v2=v2+v1*tyv(j,1)           !*tyv
        xvm3(j,i)=v1
      end do
      beta(i,1)=v2
    end do
    !$OMP END PARALLEL DO
   !*********************************************
   !print*,'beta:',beta

    !$OMP PARALLEL DO PRIVATE (v1)
    do i=1,trn
      do j=1,trn
        v1=0
        do k=1,tfn
          v1=v1+xvm3(i,k)*xvm(k,j)
        end do
        pm(i,j)=pm(i,j)-v1
      end do
    end do
    !$OMP END PARALLEL DO

    !ypy estimation
    !ypyf=MATMUL(MATMUL(transpose(tyv),pm),tyv)

    !$OMP PARALLEL DO PRIVATE (v1,v2)
    v2=0
    do i=1,trn
      v1=0
      do j=1,trn
        v1=v1+tyv(j,1)*pm(j,i)
      end do
      v2=v2+v1*tyv(i,1)
    end do
    !$OMP END PARALLEL DO

    !LKH=-0.5*(LA+LG+ypyf(1,1))
    LKH=-0.5*(LA+LG+v2)
    !print*,LA,LG,v2  !ypyf

    call cpu_time (t2_cpu)
    !print*,'LKH:',t2_cpu-t1_cpu

    if ((LKH>-99999999 .and. LKH<99999999) .and. (LKH>=LKHP)) then
      itit=1
      !print*,'ok',itit,LKHP
    else
      if (zi==1) then     !the first iteration
        sv_fact=sv_fact+0.1
        if (sv_fact > 0.99) then
          print*,'try with user-defined starting value option, -sv'
          pause
        end if
        !arbitary starting value (half of the total variance)
        vcve=0;vcva=0
        do xi=1,tn
          x1=0;x2=0
          do i=1,rn
            if (yv(i,xi).ne.-99999) then
              x1=x1+yv(i,xi)
              x2=x2+yv(i,xi)**2
            end if
          end do
          !vcve(xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))/2
          !vcva(:,xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))/(2*mn)
          vcve(xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))*(1-sv_fact)
          vcva(:,xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))*sv_fact/mn
        end do
        print*,'likelihood nan >> update with another starting values'
        goto 110
      end if

      itit=itit+1
      !LKH=-1000000000
      LKH=LKHP
      !print*,'likelihood NaN >> reduce by 0.7^',itit-1
      !print*,itit-1,'likelihood nan >> update reduced by the factor'
      print*,'likelihood nan >> update reduced by the factor'
      !print*,'not ok',itit,LKHP
      goto 111
    end if


    !PRINT '(a7,100f14.4)','LKH',LKH
    !do xi=1,tn
    !  PRINT '(a7,100f14.4)','Ve',vcve(xi,xi)
    !end do
    !print*,''
    !do wi=1,mn
    !  do xi=1,tn
    !    PRINT '(a7,100f14.4)','Va',vcva(wi,xi,xi)
    !  end do
    !  do xi=2,tn
    !    do yi=1,xi-1
    !      PRINT '(a7,100f14.4)','cov',vcva(wi,xi,yi)
    !    end do
    !  end do
    !  print*,''
    !end do
    !!pause


    ! AI matrix
    py=MATMUL(pm,tyv)       !for py

    !print*,'ai matrix'
    call cpu_time (t1_cpu)
    aim=0
    dldv=0

    nb=tn+(tn**2-tn)/2  !no block for random effects (Vg, cov)
    !Ve (diagonal)  **********************************
    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)   !total rn upto xi
      v2=0
      do i=1,rnm(xi)
        v1=0
        do j=1,rnm(xi)
          !v1=v1+py(j,1)*pm(j,i)
          !v1=v1+py(trnx-rnm(xi)+j,1)*pm(trnx-rnm(xi)+j,trnx-rnm(xi)+i)
          v1=v1+py(trnx-rnm(xi)+j,1)*resw(trnx-rnm(xi)+j)*resw(trnx-rnm(xi)+i)*pm(trnx-rnm(xi)+j,trnx-rnm(xi)+i)
        end do
        v2=v2+v1*py(trnx-rnm(xi)+i,1)
      end do 
      aim(xi,xi)=v2
    end do

    !dldv Ve ***********************************************************
    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)   !total rn upto xi
      tr1=0
      do i=1,rnm(xi)
        !tr1=tr1+pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)
        tr1=tr1+resw(trnx-rnm(xi)+i)*pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)
      end do
      tr1=-0.5*tr1
      v2=0
      do i=1,rnm(xi)
        !v2=v2+py(trnx-rnm(xi)+i,1)*py(trnx-rnm(xi)+i,1)
        v2=v2+py(trnx-rnm(xi)+i,1)*resw(trnx-rnm(xi)+i)*py(trnx-rnm(xi)+i,1)
      end do 
      dldv(xi,1)=v2*(0.5)+tr1
    end do


    !Ve (off diagonal) aim(2,1) *************************************
    trnx=rnm(1)
    do xi=2,tn
      trnx=trnx+rnm(xi)   !total rn upto xi
      trny=0
      do yi=1,(xi-1)
        trny=trny+rnm(yi)  !total rn upto yi

        v2=0
        do i=1,rnm(xi)
          v1=0
          do j=1,rnm(yi)
            !v1=v1+py(trny-rnm(yi)+j,1)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
            v1=v1+py(trny-rnm(yi)+j,1)*resw(trny-rnm(yi)+j)*resw(trnx-rnm(xi)+i)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
          end do
          v2=v2+v1*py(trnx-rnm(xi)+i,1)
        end do 
        aim(xi,yi)=v2
        !print*,v2
      end do
    end do
    !pause

    !print*,'Vg'
    !for Vg
    do k=1,mn !no. random effects

      !for Vg (doagonal) ******************************************************* 
      trnx=0
      do xi=1,tn
        trnx=trnx+rnm(xi)
        v2=0
        do i=1,rnm(xi)
          v10=0
          do j=1,rnm(xi)
            v10=v10+py(trnx-rnm(xi)+j,1)*mbin(k,yidx(xi,j),yidx(xi,i))     !py*pig
          end do
          v_tmp(trnx-rnm(xi)+i)=v10
        end do

        do i=1,trn
          v10=0
          do j=1,rnm(xi)
            v10=v10+v_tmp(trnx-rnm(xi)+j)*pm(trnx-rnm(xi)+j,i)                      !*pm
          end do
          v2_tmp(i)=v10
        end do

        do i=1,rnm(xi)
          v1=0
          do j=1,rnm(xi)
            v1=v1+v2_tmp(trnx-rnm(xi)+j)*mbin(k,yidx(xi,j),yidx(xi,i))     !*pig
          end do
          v2=v2+v1*py(trnx-rnm(xi)+i,1)                                !*py
        end do

        !aim(2+k*3-2,2+k*3-2)=v2
        !aim(2+k*3-3+xi,2+k*3-3+xi)=v2
        aim(tn+nb*k-nb+xi,tn+nb*k-nb+xi)=v2

        !for Vg x Ve (off diagonal) !******************************************
        trnv=0
        do vi=1,tn
          trnv=trnv+rnm(vi)
          v2=0
          do i=1,rnm(vi)
            !v2=v2+v2_tmp(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)    !*py
            v2=v2+v2_tmp(trnv-rnm(vi)+i)*resw(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)    !*py
          end do
          !aim(2+k*3-3+xi,vi)=v2             !with Ve1
          aim(tn+nb*k-nb+xi,vi)=v2
        end do

        !for Vg(1~[i-1]) x Vg(i)
        trnv=0
        do vi=1,xi-1
          trnv=trnv+rnm(vi)
          v2=0
          do i=1,rnm(vi)
            v1=0
            do j=1,rnm(vi)
              v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(k,yidx(vi,j),yidx(vi,i))   !*pig_tmp
            end do
            v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
          end do
          !aim(2+k*3-1,2+k*3-1-1)=v2  !wi
          !aim(2+k*3-3+xi,2+k*3-3+xi-1)=v2  !check
          aim(tn+nb*k-nb+xi,tn+nb*k-nb+vi)=v2
        end do


        !off diagonal *********************************************************
        do wi=1,k-1
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              v1=0
              do j=1,rnm(vi)
                v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(wi,yidx(vi,j),yidx(vi,i))     !*pig_tmp
              end do
              v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
            end do
            !aim(2+k*3-3+xi,2+wi*3-3+vi)=v2   !with the xi th Vg for first trait
            aim(tn+nb*k-nb+xi,tn+nb*wi-nb+vi)=v2   !with the xi th Vg for first trait
            !print*,v2
          end do

          nbc=0  !for no. block for covariance, i.e. 1: (2,1), 2: (3,1), 3:(3,2), ... 
          trnv=rnm(1)
          do vi=2,tn
            trnv=trnv+rnm(vi)   !total rn upto vi
            trnu=0
            do ui=1,(vi-1)
              nbc=nbc+1
              trnu=trnu+rnm(ui)
              v2=0
              do i=1,rnm(vi)
                v1=0
                do j=1,rnm(ui)
                  v1=v1+v2_tmp(trnu-rnm(ui)+j)*mbin(wi,yidx(ui,j),yidx(vi,i))     !*pig_tmp
                end do
                v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
              end do
              do i=1,rnm(ui)
                v1=0
                do j=1,rnm(vi)
                  v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(wi,yidx(vi,j),yidx(ui,i))     !*pig_tmp
                end do
                v2=v2+v1*py(trnu-rnm(ui)+i,1)                                !*py
              end do
              !aim(2+k*3-3+xi,2+wi*3)=v2  !check
              aim(tn+nb*k-nb+xi,tn+nb*wi-nb+tn+nbc)=v2  !check
            end do

          end do

        end do !wi

      end do  !xi
      !pause


      !dldv Vg *********************************************************
      trnx=0
      do xi=1,tn
        trnx=trnx+rnm(xi)   !total rn upto xi
        tr1=0
        do i=1,rnm(xi)
          v1=0
          do j=1,rnm(xi)
            v1=v1+pm(trnx-rnm(xi)+i,trnx-rnm(xi)+j)*mbin(k,yidx(xi,j),yidx(xi,i))
          end do
          tr1=tr1+v1
        end do
        tr1=-0.5*tr1

        v2=0
        do i=1,rnm(xi)
          v1=0
          do j=1,rnm(xi)
            v1=v1+py(trnx-rnm(xi)+j,1)*mbin(k,yidx(xi,j),yidx(xi,i))
          end do
          v2=v2+v1*py(trnx-rnm(xi)+i,1)
        end do 
        !dldv(2+k*3-3+xi,1)=v2*0.5+tr1 !************************************
        dldv(tn+nb*k-nb+xi,1)=v2*0.5+tr1 !************************************
      end do


      !for cov (diagonal) *************************************************
      nbc=0  !for no. block for covariance, i.e. 1: (2,1), 2: (3,1), 3:(3,2), ... 
      trnx=rnm(1)
      do xi=2,tn
        trnx=trnx+rnm(xi)
        trny=0
        do yi=1,xi-1
          nbc=nbc+1
          trny=trny+rnm(yi)
          v2=0
          do i=1,rnm(xi)
            v10=0
            do j=1,rnm(yi)
              v10=v10+py(trny-rnm(yi)+j,1)*mbin(k,yidx(yi,j),yidx(xi,i))     !py*pig
            end do
            v_tmp(trnx-rnm(xi)+i)=v10
          end do
          do i=1,rnm(yi)
            v10=0
            do j=1,rnm(xi)
              v10=v10+py(trnx-rnm(xi)+j,1)*mbin(k,yidx(xi,j),yidx(yi,i))     !py*pig
            end do
            v_tmp(trny-rnm(yi)+i)=v10
          end do

          !do i=1,trn     !check (may not trn)
          !  v10=0
          !  do j=1,trn
          !    v10=v10+v_tmp(j)*pm(j,i)                      !*pm
          !  end do
          !  v2_tmp(i)=v10
          !end do 

          do i=1,trn     !check (may not trn)
            v10=0
            do j=1,rnm(xi)
              v10=v10+v_tmp(trnx-rnm(xi)+j)*pm(trnx-rnm(xi)+j,i)              !*pm
            end do
            do j=1,rnm(yi)
              v10=v10+v_tmp(trny-rnm(yi)+j)*pm(trny-rnm(yi)+j,i)              !*pm
            end do
            v2_tmp(i)=v10
          end do 


          do i=1,rnm(xi)
            v1=0
            do j=1,rnm(yi)
              v1=v1+v2_tmp(trny-rnm(yi)+j)*mbin(k,yidx(yi,j),yidx(xi,i))     !*pig
            end do
            v2=v2+v1*py(trnx-rnm(xi)+i,1)                                !*py
          end do
          do i=1,rnm(yi)
            v1=0
            do j=1,rnm(xi)
              v1=v1+v2_tmp(trnx-rnm(xi)+j)*mbin(k,yidx(xi,j),yidx(yi,i))     !*pig
            end do
            v2=v2+v1*py(trny-rnm(yi)+i,1)                                !*py
          end do
          !aim(2+k*3,2+k*3)=v2     !check
          aim(tn+nb*k-nb+tn+nbc,tn+nb*k-nb+tn+nbc)=v2     !check

          !for cov x Ve (off diagonal) ************************************
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              !v2=v2+v2_tmp(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)     !*py
              v2=v2+v2_tmp(trnv-rnm(vi)+i)*resw(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)      !*py
            end do 
            !aim(2+k*3,vi)=v2          !check with Ve1
            aim(tn+nb*k-nb+tn+nbc,vi)=v2          !check with Ve1
          end do

          !for cov x Vg **************************************************
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              v1=0
              do j=1,rnm(vi)
                v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(k,yidx(vi,j),yidx(vi,i))     !*pig_tmp
              end do
              v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
            end do 
            !aim(2+k*3,2+k*3-3+vi)=v2  !with the same Vg1 for first trait
            aim(tn+nb*k-nb+tn+nbc,tn+nb*k-nb+vi)=v2  !with the same Vg1 for first trait
          end do

          !for cov x cov *************************************************
          trnv=rnm(1)
          nbc2=0    !no block for cov (vi,ui)
          do vi=2,tn
            trnv=trnv+rnm(vi)
            trnu=0
            do ui=1,(vi-1)
              nbc2=nbc2+1

              if (nbc2<nbc) then  !off diagonal within cova ************************
              trnu=trnu+rnm(ui)
              v2=0
              do i=1,rnm(vi)
                v1=0
                do j=1,rnm(ui)
                  v1=v1+v2_tmp(trnu-rnm(ui)+j)*mbin(k,yidx(ui,j),yidx(vi,i))     !*pig_tmp
                end do
                v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
              end do 
              do i=1,rnm(ui)
                v1=0
                do j=1,rnm(vi)
                  v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(k,yidx(vi,j),yidx(ui,i))     !*pig_tmp
                end do
                v2=v2+v1*py(trnu-rnm(ui)+i,1)                                !*py
              end do 
              aim(tn+nb*k-nb+tn+nbc,tn+nb*k-nb+tn+nbc2)=v2  !with the wi th cov for first trait

              end if ! *************************************************************
            end do !ui
          end do !vi


          !off diagonal *********************************************************
          do wi=1,k-1
            trnv=0
            do vi=1,tn
              trnv=trnv+rnm(vi)
              v2=0
              do i=1,rnm(vi)
                v1=0
                do j=1,rnm(vi)
                  v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(wi,yidx(vi,j),yidx(vi,i))     !*pig_tmp
                end do
                v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
              end do 
              !aim(2+k*3,2+wi*3-3+vi)=v2  !with the wi th Vg for first trait
              aim(tn+nb*k-nb+tn+nbc,tn+nb*wi-nb+vi)=v2  !with the wi th Vg for first trait
            end do

            trnv=rnm(1)
            nbc2=0    !no block for cov (vi,ui)
            do vi=2,tn
              trnv=trnv+rnm(vi)
              trnu=0
              do ui=1,(vi-1)
                nbc2=nbc2+1
                trnu=trnu+rnm(ui)
                v2=0
                do i=1,rnm(vi)
                  v1=0
                  do j=1,rnm(ui)
                    v1=v1+v2_tmp(trnu-rnm(ui)+j)*mbin(wi,yidx(ui,j),yidx(vi,i))     !*pig_tmp
                  end do
                  v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
                end do 
                do i=1,rnm(ui)
                  v1=0
                  do j=1,rnm(vi)
                    v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(wi,yidx(vi,j),yidx(ui,i))     !*pig_tmp
                  end do
                  v2=v2+v1*py(trnu-rnm(ui)+i,1)                                !*py
                end do 
                !aim(2+k*3,2+wi*3)=v2  !with the wi th cov for first trait
                aim(tn+nb*k-nb+tn+nbc,tn+nb*wi-nb+tn+nbc2)=v2  !with the wi th cov for first trait
              end do !ui
            end do !vi
          end do !wi
        end do !yi

      end do !xi **********************************************

      !dldv ***************************************************************
      nbc=0     !no block for cov (xi, yi)
      trnx=rnm(1)
      do xi=2,tn
        trnx=trnx+rnm(xi)
        trny=0
        do yi=1,xi-1
          nbc=nbc+1
          trny=trny+rnm(yi)
          tr1=0
          do i=1,rnm(xi)
            v10=0
            do j=1,rnm(yi)
              v10=v10+pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)*mbin(k,yidx(yi,j),yidx(xi,i))     !pm*pig
            end do
            tr1=tr1+v10
          end do
          do i=1,rnm(yi)
            v10=0
            do j=1,rnm(xi)
              v10=v10+pm(trny-rnm(yi)+i,trnx-rnm(xi)+j)*mbin(k,yidx(xi,j),yidx(yi,i))     !pm*pig
            end do
            tr1=tr1+v10
          end do
          tr1=-0.5*tr1

          v2=0
          do i=1,rnm(xi)
            v10=0
            do j=1,rnm(yi)
              v10=v10+py(trny-rnm(yi)+j,1)*mbin(k,yidx(yi,j),yidx(xi,i))     !py*pig
            end do
            v2=v2+v10*py(trnx-rnm(xi)+i,1)
          end do
          do i=1,rnm(yi)
            v10=0
            do j=1,rnm(xi)
              v10=v10+py(trnx-rnm(xi)+j,1)*mbin(k,yidx(xi,j),yidx(yi,i))     !py*pig
            end do
            v2=v2+v10*py(trny-rnm(yi)+i,1)
          end do
          !dldv(2+k*3,1)=v2*0.5+tr1
          dldv(tn+nb*k-nb+tn+nbc,1)=v2*0.5+tr1
        end do !yi
      end do !xi


    end do !k

    call cpu_time (t2_cpu)
    !print*,'derivatives done - time:',real(t2_cpu-t1_cpu)
    !print*,''
 
    !do i=1,(2+mn*3)
    do i=1,tn+nb*mn
      aim(i,i)=aim(i,i)/2
      do j=1,i-1
        aim(i,j)=aim(i,j)/2
        aim(j,i)=aim(i,j)
      end do
    end do

    !fixing some parameters
    wj=0
    do i=1,tn
      wj=wj+1
      if (fix_vcve(i,i) .ne. -99999) then
          aim(:,wj)=0
          aim(wj,:)=0
        aim(wj,wj)=1
      end if
    end do
    do wi=1,mn
      do xi=1,tn
        wj=wj+1
        if (fix_vcva(wi,xi,xi).ne.-99999) then
            aim(:,wj)=0
            aim(wj,:)=0
          aim(wj,wj)=1
        end if
      end do
      do xi=1,tn
        do yi=1,xi-1
          wj=wj+1
          if (fix_vcva(wi,xi,yi).ne.-99999) then
              aim(:,wj)=0
              aim(wj,:)=0
            aim(wj,wj)=1
          end if
        end do
      end do
    end do
    !*****************************

    call cholesky_inv (aim,(tn+nb*mn),x1)

    !do i=1,tn+nb*mn
    !  print*,real(aim(i,:))
    !end do
    !pause

    if (LKH-MLKH.gt.0) then
      MLKH=LKH
      !mva=va
      !mva2=va2
      !mve=ve
      !mve2=ve2
      !mcov=cov
      sdmI=aim
    end if

    !if (ABS(LKH-LKHP).lt.0.001) goto 1000
    if (ABS(LKH-LKHP).lt.conv) goto 1000

    up=MATMUL(aim,dldv)
    !do i=1,tn+nb*mn
    !  print*,real(up(i,1)),real(dldv(i,1))
    !end do
    !print*,''
    !pause

111 continue
    if (itit==1) then
      vcve_rsv=vcve
      vcva_rsv=vcva
      if (zi.ge.nit) then
        print*,'Likelihood not converged >>> may need a longer iterations'
        print*,''
        goto 1000
      end if

    else
      up=up*0.7**(itit-1)
      vcve=vcve_rsv
      vcva=vcva_rsv
    end if

    do xi=1,tn
      vcve(xi,xi)=vcve(xi,xi)+up(xi,1)
    end do
    do wi=1,mn
      do xi=1,tn
        vcva(wi,xi,xi)=vcva(wi,xi,xi)+up(tn+nb*wi-nb+xi,1)
      end do
      nbc=0
      do xi=2,tn
        do yi=1,xi-1
          nbc=nbc+1
          vcva(wi,xi,yi)=vcva(wi,xi,yi)+up(tn+nb*wi-nb+tn+nbc,1)

        end do
      end do
    end do

  END do !zi

  !va=mva;va2=mva2;ve=mve;ve2=mve2;LKH=MLKH;cov=mcov
  LKH=MLKH

1000 continue 


  open (UNIT=41,FILE=fl5,STATUS='unknown')

  do xi=1,tn
    write(41,'(a7,100f14.4)')'Ve',vcve(xi,xi),sqrt(sdmI(xi,xi))
  end do
  write(41,*)''

  do xi=1,tn
    sum_v(xi)=vcve(xi,xi)
  end do
  do wi=1,mn
    do xi=1,tn
      write(41,'(a7,100f14.4)')'Va',vcva(wi,xi,xi),sqrt(sdmI(tn+wi*nb-nb+xi,tn+wi*nb-nb+xi))
      sum_v(xi)=sum_v(xi)+vcva(wi,xi,xi)
    end do
    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        write(41,'(a7,100f14.4)')'cov',vcva(wi,xi,yi),sqrt(sdmI(tn+wi*nb-nb+tn+nbc,tn+wi*nb-nb+tn+nbc))
      end do
    end do
    write(41,*)''
  end do

  !total Var(V) *********************
  do xi=1,tn
    vvc(xi)=sdmI(xi,xi)
    do wi=1,mn
      vvc(xi)=vvc(xi)+sdmI(tn+nb*wi-nb+xi,xi)*2
      do wj=1,mn
        vvc(xi)=vvc(xi)+sdmI(tn+nb*wi-nb+xi,tn+nb*wj-nb+xi)
      end do
    end do
  end do
  !print*,vvc

  do wi=1,mn
    do xi=1,tn
      x1=vcva(wi,xi,xi)/sum_v(xi)
      x11=sdmI(tn+nb*wi-nb+xi,xi) !for the row ************

      do wj=1,mn
        x11=x11+sdmI(tn+nb*wi-nb+xi,tn+nb*wj-nb+xi)
      end do
      !print*,x11

      !SE **********************************
      x2=x1**2*(sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+xi)/vcva(wi,xi,xi)**2+vvc(xi)/sum_v(xi)**2-2*x11/(vcva(wi,xi,xi)*sum_v(xi)))
      write(41,'(a7,100f14.4)')'h2',x1,sqrt(x2)

    end do
    write(41,*)''

    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        z1=vcva(wi,xi,yi)/sqrt(vcva(wi,xi,xi)*vcva(wi,yi,yi))
        !var(cor) *************************************
        z2=sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+xi)/(4*vcva(wi,xi,xi)**2)
        z2=z2+sdmI(tn+nb*wi-nb+yi,tn+nb*wi-nb+yi)/(4*vcva(wi,yi,yi)**2)
        z2=z2+sdmI(tn+nb*wi-nb+tn+nbc,tn+nb*wi-nb+tn+nbc)/(vcva(wi,xi,yi)**2)
        z2=z2+2*sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+yi)/(4*vcva(wi,xi,xi)*vcva(wi,yi,yi))
        z2=z2-2*sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+tn+nbc)/(2*vcva(wi,xi,xi)*vcva(wi,xi,yi))
        z2=z2-2*sdmI(tn+nb*wi-nb+tn+nbc,tn+nb*wi-nb+yi)/(2*vcva(wi,xi,yi)*vcva(wi,yi,yi))
        z2=z2*(z1**2)
        write(41,'(a7,100f14.4)')'cor',z1,sqrt(z2)
      end do
    end do
    write(41,*)''
  end do

  write(41,*)''
  write(41,'(a7,100f14.4)')'LKH',LKH

  write(41,*)''
  do i=1,(tn+nb*mn)
    !write(41,*)real(sdmI(i,1:i))
    !write(41,'(100f12.5)')real(sdmI(i,1:i))
    write(41,'(100g18.8)')(sdmI(i,1:i))
  end do

  close(41)

!CDF functions *****************************************************
      !p = 1D+00 - k !0.011135D+00
      !q = 1.0D+00 - p
      p = huge ( p )
      q = huge ( q )
      x = beta(tfn,1)**2/xvmx(tfn,tfn)
      !x=-1.28155 !2 3.0D+00
      !mean = 0.0D+00 ! 5.0D+00
      sd = 1.0D+00 !0.875D+00    !degree of greedom for cdfchi
which=1  !Calculate P and Q from X, MEAN and SD;
!which=2  !Calculate X from P, Q, MEAN and SD;
!which=3  !Calculate MEAN from P, Q, X and SD;
!which=4  !Calculate SD from P, Q, X and MEAN;

call cdfchi ( which, p, q, x, sd, status, bound )


  !do xi=1,tfn
  !  write(47,'(i3,100f24.16)')xi,beta(xi,1),sqrt(xvmx(xi,xi))
  !end do
  write(47,'(a3,a17,a3,i8,2f14.5,f12.3,e14.5)')cno(gwasi),trim(markt(gwasi)),mafreq1(gwasi,1),nmiss(gwasi),real(beta(tfn,1)),real(sqrt(xvmx(tfn,tfn))),real(beta(tfn,1)**2/xvmx(tfn,tfn)),real(q)




end do !gwasi


close(47)


end subroutine

!end module
