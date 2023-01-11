
subroutine aireml_m_h_sblup (mbin,yidx,yv,rn,rnm,trn,pedn,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,fl24,xmm,nit,conv,wt_res,wmj,fl111,fl111r)        

!***********************************************************************
!BLUP and SNP BLUP  
!S. Hong Lee (2020)

!mbin : matrices bin
!obs  : observed ID
!yv   : phenotypes
!rn   : no. phenotypes
!pedn : no. pedigree (ID) should be = rn
!fn   : no. fixed effects
!mn   : no. matrices (random effects)
!wm   : a pedn x wmi matrix having variables (wmi is # variables)
!***********************************************************************

implicit none

INTEGER::n,nped,ix,iz,tn,nb,nbc,nbc2,xi,yi,ri,zi,vi,ui,wi,wmj,wmi,wj,nit
integer::rn,pedn,mn,trn,trnx,trny,tfn,trnv,trnu,itit,nm,im
integer::yidx(tn,rn),rnm(tn),fnm(tn),eoro
integer::rn1,rn2,yidx1(rn),yidx2(rn)
real::mbin(mn,pedn,pedn),yv(rn,tn),sv_fact,sclf
real,allocatable::wm(:,:),mafreq2(:,:),mpos(:),htz(:,:)
double precision::resw(trn),resw2(pedn,tn),tyv(trn,1)

character(len=1),allocatable::mafreq1(:,:)
character(len=3),allocatable::cno(:)
character(len=100),allocatable::pedor(:,:)
character(len=50),allocatable::markt(:)

double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),conv
double precision::vcva(mn,tn,tn),vcve(tn,tn),vvc(tn)
double precision::fix_vcva(mn,tn,tn),fix_vcve(tn,tn)
double precision::vcva_rsv(mn,tn,tn),vcve_rsv(tn,tn)
double precision::blup_ebv(pedn,tn,mn),blup_py(pedn,tn,mn),beta(tfn,1) !,cor
double precision::blup_r(pedn,tn,mn),sblup_v,sblup_pev
double precision,allocatable::sblup_ebv(:,:),sblup_r(:,:)
double precision,allocatable::tmp_v(:,:)

integer::i,j,m,k,l,io,plink_file
double precision::sum_h2,sum_v(tn)

double precision::LC,LA,LG,ypy,x1,x2,x3,x10,x11,v1,v2,v10,v11
double precision::y1,y2,y3,y10,y11,z1,z2,z3,z10,z11
double precision::LKHP,MLKH,mva(mn),mva2(mn),mve,mve2,mcov(mn)
double precision::h2(mn),h2n(mn),tr1,tr2(mn),tr3

double precision::xm(trn,tfn),xmm(tn,pedn,10000)

double precision::xvmx(tfn,tfn),xvm(tfn,trn) 
double precision::xvm2(tfn,trn),xvm3(trn,tfn) 
double precision::pm(trn,trn),py(trn,1),xb(trn,1)

!double precision::v_tmp(trn),v2_tmp(trn)

!double precision::aim(tn+mn*(tn+(tn**2-tn)/2),tn+mn*(tn+(tn**2-tn)/2))
!double precision::sdmI(tn+mn*(tn+(tn**2-tn)/2),tn+mn*(tn+(tn**2-tn)/2))
!double precision::up(tn+mn*(tn+(tn**2-tn)/2),1),dldv(tn+mn*(tn+(tn**2-tn)/2),1)

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
character(:),allocatable::string2
character(len=64)::fl4,fl24,fl5,fl11,fl11r,fl111,fl111r,wt_res,cdum_t(tn)
character(len=64)::filnam
logical::file_exist

!CDF
real ( kind = 8 ) bound,mean,p,q,sd,x
integer ( kind = 4 ) status,which

!sequential storing
integer*1,allocatable::sqgenom(:,:),sex(:)
integer*1 b1,igen(0:1,0:1),w1,w2

xm=0;trn=0;tfn=0
do xi=1,tn
  tfn=tfn+fnm(xi)
  do i=1,rn
    if (yv(i,xi).ne.-99999) then
      trn=trn+1
      tyv(trn,1)=yv(i,xi)
      do j=1,fnm(xi)
        xm(trn,tfn-fnm(xi)+j)=xmm(xi,i,j)
      end do
    end if
  end do
end do


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

110 continue

  itit=1
  ! Start iteration
  LKH=-1000000000
  MLKH=-1000000000


  ! V matrix    V= ZAZ'Va+Z2GZ2'Vpe+IVe*************************
  !do zi=1,nit
  do zi=1,1
    LKHP=LKH
  !LKH and AIREML (without sparce technique)
  ! V matrix    V = ZAZ'Va + IVe *************************

  call cpu_time (t2_cpu)

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
          end do
          if (xi.ne.yi) then
            pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)
          end if
        end do
        if (xi==yi) then
          pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)+vcve(xi,xi)*resw(trnx-rnm(xi)+i)
        end if
      end do
      !$OMP END PARALLEL DO
    end do
  end do
  call cpu_time (t2_cpu)
  call cholesky_inv (pm,trn,LA)
  call cpu_time (t1_cpu)
  print*,'V inverse done - time:',real(t1_cpu-t2_cpu) !,LA

  ! P = (VI)-(VI X (X'VI X)I' X' VI)******************************
  !xvm=MATMUL(transpose(xm),pm)
  call dgemm ('T','N',tfn,trn,trn,1.0D0,xm,trn,pm,trn,0.0D0,xvm,tfn)
  call dgemm ('N','N',tfn,tfn,trn,1.0D0,xvm,tfn,xm,trn,0.0D0,xvmx,tfn)

  call cholesky_inv (xvmx,tfn,LG)

  call dgemm ('T','N',trn,tfn,tfn,1.0D0,xvm,tfn,xvmx,tfn,0.0D0,xvm3,trn)
  call dgemm ('T','N',tfn,1,trn,1.0D0,xvm3,trn,tyv,trn,0.0D0,beta,tfn) ! tyv is real?

  call dgemm('N','N',trn, trn, tfn, -1.0D0, xvm3, trn, xvm, tfn, 1.0D0,pm,trn)

    !ypy estimation
    !ypyf=MATMUL(MATMUL(transpose(tyv),pm),tyv)
    allocate(tmp_v(trn,1))
    !$OMP PARALLEL DO PRIVATE (i) 
    do i=1,trn
      tmp_v(i,1)=dot_product(tyv(:,1),pm(:,i))
    end do
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO PRIVATE (i) 
    do i=1,trn
      v2=dot_product(tyv(:,1),tmp_v(:,1))
    end do
    !$OMP END PARALLEL DO
    deallocate(tmp_v)

    LKH=-0.5*(LA+LG+v2)
    !print*,LA,LG,v2  !ypyf

    call cpu_time (t2_cpu)
    !print*,'LKH:',t2_cpu-t1_cpu

    if ((LKH>-99999999 .and. LKH<99999999) .and. (LKH>=LKHP)) then
      itit=1
    else
      print*,'NaN LKH, something wrong >>> check starting value and input'
      pause
    end if

    PRINT '(a7,100f14.4)','LKH',LKH
    do xi=1,tn
      PRINT '(a7,100f14.4)','Ve',vcve(xi,xi)
    end do
    print*,''
    do wi=1,mn
      do xi=1,tn
        PRINT '(a7,100f14.4)','Va',vcva(wi,xi,xi)
      end do
      do xi=2,tn
        do yi=1,xi-1
          PRINT '(a7,100f14.4)','cov',vcva(wi,xi,yi)
        end do
      end do
      print*,''
    end do
    !pause

    ! AI matrix
    py=MATMUL(pm,tyv)       !for py


  END do !zi

if (fl11.ne.'null') then

  print*,'for individual BLUP solution with the starting value ************'

  !make it triangle -> full ****************
  do xi=1,tn
    do yi=1,xi
      do k=1,mn
        vcva(k,yi,xi)=vcva(k,xi,yi)
        !print*,vcva
      end do
    end do
  end do !***********************************

  trnx=0
  do xi=1,tn
    trnx=trnx+rnm(xi)
    do zi=1,pedn
      do k=1,mn
        v1=0;v2=0
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          do j=1,rnm(yi)
            v1=v1+mbin(k,yidx(yi,j),zi)*vcva(k,xi,yi)*py(trny-rnm(yi)+j,1)
            if (yidx(yi,j)==zi) then
              v2=v2+vcva(k,xi,yi)*py(trny-rnm(yi)+j,1)
            end if
          end do
        end do
        blup_ebv(zi,xi,k)=v1
        blup_py(zi,xi,k)=v2
      end do !k
    end do !zi
  end do !xi

  open (unit=45,file=trim(fl11)//".py",status='unknown')
  do xi=1,tn
    do i=1,pedn
      write(45,'(i3,100f24.16)')xi,blup_py(i,xi,:)
    end do
  end do
  close(45)

  open (unit=47,file=trim(fl11)//".fsl",status='unknown')
  write(47,'(a3,2a14,a12,a14)')'#','BETA','SE','CHI','P'

  do xi=1,tfn
    !CDF functions *****************************************************
    p = huge ( p )
    q = huge ( q )
    x = beta(xi,1)**2/xvmx(xi,xi)
    sd = 1.0D+00 !degree of greedom for cdfchi
    which=1  !Calculate P and Q from X, MEAN and SD;
    call cdfchi ( which, p, q, x, sd, status, bound )
    write(47,'(i3,2f14.5,f12.3,e14.5)')xi,beta(xi,1),sqrt(xvmx(xi,xi)),x,q
  end do
  close(47)

  open (unit=43,file=fl11,status='unknown')
  !write(43,'(100f24.16)')beta
  write(43,'(a5,a26)')'trait','EBVs ...'
  do xi=1,tn
    do i=1,pedn
      write(43,'(i3,100f24.16)')xi,blup_ebv(i,xi,:)
    end do
  end do
  close(43)

end if

if (fl11r.ne.'null') then
  print*,'for BLUP solution with the starting value ***********'
  !make it triangle -> full ****************
  do xi=1,tn
    do yi=1,xi
      do k=1,mn
        vcva(k,yi,xi)=vcva(k,xi,yi)
        !print*,vcva
      end do
    end do
  end do !***********************************

  trnx=0
  do xi=1,tn
    trnx=trnx+rnm(xi)
    do zi=1,pedn
      do k=1,mn
        v1=0;v2=0
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          do j=1,rnm(yi)
            v1=v1+mbin(k,yidx(yi,j),zi)*vcva(k,xi,yi)*py(trny-rnm(yi)+j,1)
            if (yidx(yi,j)==zi) then
              v2=v2+vcva(k,xi,yi)*py(trny-rnm(yi)+j,1)
            end if
          end do
        end do
        blup_ebv(zi,xi,k)=v1
        blup_py(zi,xi,k)=v2
      end do !k
    end do !zi
  end do !xi

  !to get reliability for BLUP, i.e. GPA
  !$OMP PARALLEL DO PRIVATE(vi,zi,k,v1,v2)
  do vi=1,tn
    do zi=1,pedn
      do k=1,mn
        v2=0
        trnx=0
        do xi=1,tn
          trnx=trnx+rnm(xi)
          do i=1,rnm(xi)
            v1=0
            trny=0
            do yi=1,tn
              trny=trny+rnm(yi)
              do j=1,rnm(yi)
                v1=v1+mbin(k,yidx(yi,j),zi)*vcva(k,vi,yi)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
              end do
            end do
            v2=v2+v1*mbin(k,yidx(xi,i),zi)
          end do
        end do
        blup_r(zi,vi,k)=v2
      end do
    end do
  end do
  !$OMP END PARALLEL DO 

  open (unit=45,file=trim(fl11r)//".py",status='unknown')
  do xi=1,tn
    do i=1,pedn
      write(45,'(i3,100f24.16)')xi,blup_py(i,xi,:)
    end do
  end do
  close(45)

  open (unit=47,file=trim(fl11r)//".fsl",status='unknown')
  write(47,'(a3,2a14,a12,a14)')'#','BETA','SE','CHI','P'
  do xi=1,tfn
    !CDF functions *****************************************************
    p = huge ( p )
    q = huge ( q )
    x = beta(xi,1)**2/xvmx(xi,xi)
    sd = 1.0D+00 !degree of greedom for cdfchi
    which=1  !Calculate P and Q from X, MEAN and SD;
    call cdfchi ( which, p, q, x, sd, status, bound )
    write(47,'(i3,2f14.5,f12.3,e14.5)')xi,beta(xi,1),sqrt(xvmx(xi,xi)),x,q
  end do
  close(47)

  open (unit=43,file=fl11r,status='unknown')
  !write(43,'(100f24.16)')beta
  write(43,'(a5,a26)')'trait','EBVs ...'
  do xi=1,tn
    do i=1,pedn
      write(43,'(i3,100f24.16)')xi,blup_ebv(i,xi,:)
    end do
  end do
  close(43)

  open (unit=46,file=trim(fl11r)//".r2",status='unknown')
  do xi=1,tn
    do i=1,pedn
      !write(46,'(i3,100f24.16)')xi,blup_r(i,xi,:)**2
      write(46,'(i3,100f24.16)')xi,blup_r(i,xi,:)
    end do
  end do
  close(46)

end if

!SNP BLUP
!W mat or PLINK files *****************************************
if (fl111.ne.'null') then
  filnam=fl111
elseif (fl111r.ne.'null') then
  filnam=fl111r
end if

  inquire(file=filnam,exist=file_exist)
  if (file_exist) then
    plink_file=0
  else
    inquire(file=trim(filnam)//".bim",exist=file_exist)
    if (file_exist) then
      inquire(file=trim(filnam)//".bed",exist=file_exist)
      if (file_exist) then
        inquire(file=trim(filnam)//".fam",exist=file_exist)
        if (file_exist) then
          plink_file=1
        else
          print*,'plink fam file missing'
          pause
        end if
      else
        print*,'plink bed file missing'
        pause
      end if
    else
      print*,'plink bed file missing'
      pause
    end if
  end if
  !********************************************************

if (plink_file==1) then
  nm=0
  open (unit=36,file=trim(filnam)//".bim",status='old')
  do
    nm=nm+1
    read(36,*,iostat=io)cdum
    if (io.ne.0) exit
  end do
  close(36)
  nm=nm-1
  print*,'no. marker: ',nm
  wmi=nm

  allocate(mpos(nm),cno(nm),mafreq1(nm,2),mafreq2(nm,1))
  allocate(pedor(pedn,4),markt(nm),htz(nm,1))
  allocate(sqgenom(int(pedn/4)+1,nm),sex(pedn))

  open (unit=38,file=trim(filnam)//'.bed',status='old',access='stream',form='unformatted')
  open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
  open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file

  do i=1,nm
    read(36,*)cno(i),markt(i),cdum,mpos(i),mafreq1(i,1),mafreq1(i,2)
  end do

!for fam file (ped) ***************************************************
do i=1,pedn
  read(35,*)pedor(i,1:4),sex(i),cdum
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
  do i=1,int((pedn-1)/4)+1    !no. ID / 4
    read(38)sqgenom(i,im)
  end do
  !print*,igenom(i,1:20)
  !pause

end do
write (*,'(a22)')'100% read plink files'    !closing progress report

close(35)
close(38)
close(36)

!check duplications in ped file **********************************
do i=1,pedn
  do j=i+1,pedn
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

eoro=2    !check
do zi=1,1
  PRINT*,'allele frequency' ! for pop',zi
  !$OMP PARALLEL DO PRIVATE (yi,xi,w1) 
  do im=1,nm
    x1=0;x2=0;x3=0
    do i=1,pedn
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
    !print*,htz(im,zi)
  end do
  !$OMP END PARALLEL DO
end do

sclf=-1       !check
allocate(wm(pedn,wmi))
  !$OMP PARALLEL DO PRIVATE (yi,xi,w1,x3) 
  do im=1,nm
    do i=1,pedn
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
  !$OMP END PARALLEL DO
  !print*,wm(1,1:5)
  !print*,wm(2,1:5)
  !print*,wm(3,1:5)

else
    allocate(character(100000000)::string2)
    open (unit=71,file=filnam,status='old')
    read(71,'(a)')string2
    close(71)
    do i=1,100000000
      read(string2,*,iostat=io)(cdum,j=1,i)
      if (io.ne.0) exit
    end do
    wmi=i-1
    print*,'# variables:',wmi
    deallocate(string2)
  allocate(wm(pedn,wmi))

    open (unit=71,file=filnam,status='old')
    do i=1,pedn
      read(71,*)(wm(i,j),j=1,wmi)
    end do
    close(71)

    allocate(mafreq1(wmi,2))
    mafreq1="0"
end if   !plink_file == 1

allocate(sblup_ebv(wmi,tn))

  !make it triangle -> full ****************
  do xi=1,tn
    do yi=1,xi
      do k=1,mn
        vcva(k,yi,xi)=vcva(k,xi,yi)
        !print*,vcva
      end do
    end do
  end do !***********************************

  !for ith specific component (i.e. -snpblup i (i=wmj))
  trnx=0
  do xi=1,tn
    trnx=trnx+rnm(xi)
    !do zi=1,pedn
    do zi=1,wmi
      !do k=1,mn
      do k=wmj,wmj
        v1=0;v2=0
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          do j=1,rnm(yi)
            !v1=v1+mbin(k,yidx(yi,j),zi)*vcva(k,xi,yi)*py(trny-rnm(yi)+j,1)
            v1=v1+wm(yidx(yi,j),zi)*(vcva(k,xi,yi)/wmi)*py(trny-rnm(yi)+j,1)
          end do
        end do
        sblup_ebv(zi,xi)=v1
      end do !k
    end do !zi
  end do !xi

if (fl111.ne.'null') then
  print*,'for SNP BLUP solution with starting values *******************'

  open (unit=48,file=trim(fl111)//".fsl",status='unknown')
  write(48,'(a3,2a14,a12,a14)')'#','BETA','SE','CHI','P'

  do xi=1,tfn
    !CDF functions *****************************************************
    p = huge ( p )
    q = huge ( q )
    x = beta(xi,1)**2/xvmx(xi,xi)
    sd = 1.0D+00 !degree of greedom for cdfchi
    which=1  !Calculate P and Q from X, MEAN and SD;
    call cdfchi ( which, p, q, x, sd, status, bound )
    write(48,'(i3,2f14.5,f12.3,e14.5)')xi,beta(xi,1),sqrt(xvmx(xi,xi)),x,q
  end do
  close(48)

  open (unit=49,file=trim(fl111)//".sumstat",status='unknown')
  write(49,'(a5,a7,a13)')'trait','Ref_A','SNP EBVs'
  do xi=1,tn
    !do i=1,pedn
    do i=1,wmi
      !write(49,'(i3,1f14.5)')xi,sblup_ebv(i,xi)
      write(49,'(i5,a5,1f14.5)')xi,mafreq1(i,1),sblup_ebv(i,xi)
    end do
  end do
  close(49)

elseif (fl111r.ne.'null') then
  print*,'for SNP BLUP solution,r2,SE,p-value with starting values ********'

  allocate(sblup_r(wmi,tn))

  !to get reliability for SNP BLUP, i.e. W'(va/mv)PW see se3.r in mac/stuff/erm
  !$OMP PARALLEL DO PRIVATE(vi,zi,k,v1,v2)
  do vi=1,tn
    !do zi=1,pedn
    do zi=1,wmi
      !do k=1,mn
      do k=wmj,wmj
        v2=0
        trnx=0
        do xi=1,tn
          trnx=trnx+rnm(xi)
          do i=1,rnm(xi)
            v1=0
            trny=0
            do yi=1,tn
              trny=trny+rnm(yi)
              do j=1,rnm(yi)
                !v1=v1+mbin(k,yidx(yi,j),zi)*vcva(k,vi,yi)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
                v1=v1+wm(yidx(yi,j),zi)*(vcva(k,vi,yi)/wmi)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
              end do
            end do
            !v2=v2+v1*mbin(k,yidx(xi,i),zi)
            v2=v2+v1*wm(yidx(xi,i),zi)
          end do
        end do
        sblup_r(zi,vi)=v2
      end do
    end do
  end do
  !$OMP END PARALLEL DO 

  open (unit=48,file=trim(fl111r)//".fsl",status='unknown')
  write(48,'(a3,2a14,a12,a14)')'#','BETA','SE','CHI','P'

  do xi=1,tfn
    !CDF functions *****************************************************
    p = huge ( p )
    q = huge ( q )
    x = beta(xi,1)**2/xvmx(xi,xi)
    sd = 1.0D+00 !degree of greedom for cdfchi
    which=1  !Calculate P and Q from X, MEAN and SD;
    call cdfchi ( which, p, q, x, sd, status, bound )
    write(48,'(i3,2f14.5,f12.3,e14.5)')xi,beta(xi,1),sqrt(xvmx(xi,xi)),x,q
  end do
  close(48)

  open (unit=49,file=trim(fl111r)//".sumstat",status='unknown')
  !write(49,'(a5,4a13,2a13)')'trait','SNP_EBVs','r2','SE','PEV','Chi','P'
  write(49,'(a5,a7,4a13,2a13)')'trait','Ref_A','SNP_EBVs','r2','SE','PEV','Chi','P'
  do xi=1,tn
    !do i=1,pedn
    do i=1,wmi
      sblup_v=sblup_r(i,xi)*(vcva(wmj,xi,xi)/wmi)
      sblup_pev=(1-sblup_r(i,xi))*(vcva(wmj,xi,xi)/wmi)
      x = sblup_ebv(i,xi)**2/sblup_v
      call cdfchi ( which, p, q, x, sd, status, bound )
      !write(49,'(i3,100f24.16)')xi,sblup_ebv(i,xi),sblup_r(i,xi),sqrt(sblup_v),sblup_pev,x,q
      !write(49,'(i3,4f14.5,f12.3,e14.5)')xi,real(sblup_ebv(i,xi)),real(sblup_r(i,xi)),real(sqrt(sblup_v)),real(sblup_pev),real(x),real(q)
      write(49,'(i5,a5,4f14.5,f12.3,e14.5)')xi,mafreq1(i,1),real(sblup_ebv(i,xi)),real(sblup_r(i,xi)),real(sqrt(sblup_v)),real(sblup_pev),real(x),real(q)
    end do
  end do
  close(49)

  !deallocate(wm,sblup_ebv,sblup_r)
end if



end subroutine

