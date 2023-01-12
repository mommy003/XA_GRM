
subroutine aireml_rxe (yidx,yv,rn,rnm,trn,pedn,ttn,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,fl12,fl24,xmm,nit,conv,thn)        

!***********************************************************************
!aireml_rxe: RxE estimation 
!S. Hong Lee (2019)

!obs  : observed ID
!yv   : phenotypes
!rn   : no. phenotypes
!pedn : no. pedigree (ID) should be = rn
!fn   : no. fixed effects 
!mn   : no. matrices (random effects)
!***********************************************************************

implicit none

INTEGER::n,nped,ix,iz,tn,nb,nb2,nb3,nbc,nbc2,xi,yi,yi2,ri,zi,vi,ui,vi2,ui2,wi,wj
integer::ttn,m_tnk(ttn,mn),sub_tn(ttn)    !no. polynomial components
integer::rn,pedn,mn,trn,trnx,trny,tfn,trnv,trnu,itit,nit
integer::yidx(tn,rn),rnm(tn),fnm(tn)
integer::rn1,rn2,ll,mm,mm2,thn
!real::mbin(mn,pedn,pedn),yv(rn,tn),tyv(trn,1)
real::yv(rn,tn),tyv(trn,1),res(trn,1),d_beta(tfn,1),d_res(trn,1)

double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),conv
double precision::vcva(mn,trn),pig_vcva(mn,trn),vcve(tn,tn),fix_vcve(tn,tn)
double precision::vcve_rsv(tn,tn)
double precision::blup_ebv(pedn,100,mn),beta(tfn,1) ! check total sum k (pol oreder) < 100
double precision::blup_r(pedn,100,mn) ! total sum k (pol oreder) < 100

integer::i,j,m,k,k2,l,io
double precision::sum_h2,sum_v(tn)

double precision::LC,LA,LG,ypy,x1,x2,x3,x10,x11,v1,v2,v3,v10,v11
double precision::y1,y2,y3,y10,y11,z1,z2,z3,z10,z11
double precision::LKHP,MLKH,mva(mn),mva2(mn),mve,mve2,mcov(mn)
double precision::h2(mn),h2n(mn),tr1,tr2(mn),tr3

double precision::xm(trn,tfn),xmm(tn,pedn,1000)

double precision::xvmx(tfn,tfn),xvm(tfn,trn) 
double precision::xvm2(tfn,trn),xvm3(trn,tfn),tr_vmi
!double precision::pm(trn,trn),pm2(trn),py(trn,1),xb(trn,1)
double precision::vmi(trn),pm2(trn),py(trn,1),xb(trn,1),llik(trn),llik2(trn)
double precision::pig(trn),pig_tmp(trn),dpy(trn,1),xvmdpy(tfn,1),pdpy(trn,1)

double precision::v_tmp(trn),v2_tmp(trn)

!time measure **********************************************************
INTEGER::now(8)
CHARACTER*8 date
CHARACTER*20 time,zone
double precision:: timef, t1_real, t2_real, elapsed_real_secs
double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
!***********************************************************************

character(len=24)::cdum
character(len=64)::fl4,fl5,fl11,fl11r,fl12,fl24
logical::file_exist

double precision::m_tkm(mn,ttn,20,20),m_tkm_rsv(mn,ttn,20,20),fix_m_tkm(mn,ttn,20,20)
double precision::ctkm(mn,ttn,ttn,20,20),ctkm_rsv(mn,ttn,ttn,20,20),fix_ctkm(mn,ttn,ttn,20,20)
double precision,allocatable::km(:,:),km2(:,:),phi(:,:),km2_2(:,:),phi_2(:,:)
double precision,allocatable::phi2(:,:),phi_1(:,:),tmp_v(:,:)
!double precision::wk(tn),xk(tn)
double precision::wk(pedn),m_wk(ttn,pedn)
double precision,allocatable::tmp_wk(:),tmp_wk2(:),tmp_wk_2(:),tmp_wk_1(:)
double precision,allocatable::km_phi(:,:),km_phi2(:,:),gr_km_phi(:,:)

double precision,allocatable::aim(:,:),infm(:,:)
double precision,allocatable::sdmI(:,:)
double precision,allocatable::up(:,:),dldv(:,:),dldv_sc(:,:)

integer::idum,idum2,ii,jj,kk,oo,ii2,jj2,xi2,i2,j2

!CDF
real ( kind = 8 ) bound,mean,p,q,sd,x
integer ( kind = 4 ) status,which

open (UNIT=46,FILE=fl12,STATUS='old')
read(46,*)(sub_tn(i),i=1,ttn)
do i=1,ttn
  if (sub_tn(i).ne.rnm(i)) then
    print*,'# E variable and # phenotype not matched',i,rnm(i),sub_tn(i)
    pause
  end if
end do

do i=1,ttn
  read(46,*)m_tnk(i,:)
  !read(46,*)m_wk(i,1:sub_tn(i))
  !do j=1,sub_tn(i)
  k=0
  do j=1,rn
    read(46,*,iostat=io)cdum
    if (io.ne.0) then
      print*,'# E variable not matched with # non-missing phenotype',j
      pause
    end if
    if (yv(j,i).ne.-99999 .and. cdum=='NA') then
      print*,'missing E variable for non-missing phenotype',j
      pause
    end if
    if (yv(j,i).ne.-99999) then
      k=k+1
      read(cdum,*)m_wk(i,k)
    end if
  end do
  if (k.ne.sub_tn(i)) then
    print*,'E variable and phenotype not matched',i,k,sub_tn(i)
    pause
  end if
end do
close(46)

idum=tn
do ll=1,mn
  do kk=1,ttn
    idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
  end do
  do kk=2,ttn
    do mm=1,kk-1
      idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
    end do
  end do
end do

allocate(aim(idum,idum),infm(idum,idum),sdmI(idum,idum),up(idum,1),dldv(idum,1),dldv_sc(rn,idum))

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

!xm=1  !check
!tfn=1  !check

print*,'*** number of records used ***'
do i=1,tn
  print*,'site',i,':',rnm(i)
end do
print*,''


  ! y = Xb + Zu + e ---- fitting effects without QTL
  !arbitary starting value (half of the total variance)
  vcve=0;m_tkm=0;ctkm=0   !;km=0
  do xi=1,tn
    x1=0;x2=0
    do i=1,rn
      if (yv(i,xi).ne.-99999) then
        x1=x1+yv(i,xi)
        x2=x2+yv(i,xi)**2
      end if
    end do
    vcve(xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))/2
    !vcva(:,xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))/(2*mn)
  end do
  
  inquire(file=fl4,exist=file_exist)
  if (file_exist) then
    open (UNIT=44,FILE=fl4,STATUS='old')
    do i=1,tn
      read(44,*,iostat=io)cdum,vcve(i,i)
      if (io.ne.0) exit
    end do
    do wi=1,mn
      do ii=1,ttn
        allocate(km(m_tnk(ii,wi),m_tnk(ii,wi)))

        do xi=1,m_tnk(ii,wi)
          read(44,*,iostat=io)cdum,km(xi,xi)
          if (io.ne.0) exit
        end do
        do xi=1,m_tnk(ii,wi)
          do yi=1,xi-1
            read(44,*,iostat=io)cdum,km(xi,yi)
            km(yi,xi)=km(xi,yi)
            if (io.ne.0) exit
          end do
        end do
        m_tkm(wi,ii,1:m_tnk(ii,wi),1:m_tnk(ii,wi))=km
        deallocate(km)
      end do
      do ii=2,ttn
        do jj=1,ii-1
          allocate(km(m_tnk(ii,wi),m_tnk(jj,wi)))
          do xi=1,m_tnk(ii,wi)
            do yi=1,m_tnk(jj,wi)
              read(44,*,iostat=io)cdum,km(xi,yi)
              if (io.ne.0) exit
            end do
          end do
          ctkm(wi,ii,jj,1:m_tnk(ii,wi),1:m_tnk(jj,wi))=km
          deallocate(km)
        end do
      end do

    end do !wi
    close(44)
  end if

  inquire(file=fl24,exist=file_exist)
  if (file_exist) then
    open (UNIT=51,FILE=fl24,STATUS='old')
    do i=1,tn
      read(51,*,iostat=io)cdum,fix_vcve(i,i)
      if (io.ne.0) exit
    end do
    do wi=1,mn
      do ii=1,ttn
        allocate(km(m_tnk(ii,wi),m_tnk(ii,wi)))

        do xi=1,m_tnk(ii,wi)
          read(51,*,iostat=io)cdum,km(xi,xi)
          if (io.ne.0) exit
        end do
        do xi=1,m_tnk(ii,wi)
          do yi=1,xi-1
            read(51,*,iostat=io)cdum,km(xi,yi)
            km(yi,xi)=km(xi,yi)
            if (io.ne.0) exit
          end do
        end do
        fix_m_tkm(wi,ii,1:m_tnk(ii,wi),1:m_tnk(ii,wi))=km
        deallocate(km)
      end do
      do ii=2,ttn
        do jj=1,ii-1
          allocate(km(m_tnk(ii,wi),m_tnk(jj,wi)))
          do xi=1,m_tnk(ii,wi)
            do yi=1,m_tnk(jj,wi)
              read(51,*,iostat=io)cdum,km(xi,yi)
              if (io.ne.0) exit
            end do
          end do
          fix_ctkm(wi,ii,jj,1:m_tnk(ii,wi),1:m_tnk(jj,wi))=km
          deallocate(km)
        end do
      end do

    end do !wi
    close(51)
  else
    fix_vcve=-99999
    fix_m_tkm=-99999
    fix_ctkm=-99999
  end if

  itit=1   !iteration in the iteration

  ! Start iteration
  LKH=-1000000000
  MLKH=-1000000000

  print*,'iteration start'
  do zi=1,10000000
  call cpu_time (t2_cpu)
    LKHP=LKH
  !LKH and AIREML (without sparce technique)
  !polynomial components not same for all multiple random effects

    !fixing some parameters
    do i=1,tn
      if (fix_vcve(i,i) .ne. -99999) vcve(i,i)=fix_vcve(i,i)
    end do
    do wi=1,mn
      do ii=1,ttn
        do xi=1,m_tnk(ii,wi)
          do yi=1,m_tnk(ii,wi)
            if (fix_m_tkm(wi,ii,xi,yi).ne.-99999) then
              m_tkm(wi,ii,xi,yi)=fix_m_tkm(wi,ii,xi,yi)
            end if
          end do
        end do
      end do
      do ii=2,ttn
        do jj=1,ii-1
          do xi=1,m_tnk(ii,wi)
            do yi=1,m_tnk(jj,wi)
              if (fix_ctkm(wi,ii,jj,xi,yi).ne.-99999) then
                ctkm(wi,ii,jj,xi,yi)=fix_ctkm(wi,ii,jj,xi,yi)
              end if
            end do
          end do
        end do
      end do
    end do


  do wi=1,mn
    do ii=1,ttn
      allocate(phi(rnm(ii),m_tnk(ii,wi)),km(m_tnk(ii,wi),m_tnk(ii,wi)))
      allocate(tmp_wk(rnm(ii)))

      tmp_wk=m_wk(ii,1:rnm(ii))
      call pol (m_tnk(ii,wi),tmp_wk,rnm(ii),phi)

      km=m_tkm(wi,ii,1:m_tnk(ii,wi),1:m_tnk(ii,wi))

      mm=sum(rnm(1:ii))-rnm(ii)
      !$OMP PARALLEL DO PRIVATE(jj)
      do jj=1,rnm(ii)
        v1=0
        do i=1,m_tnk(ii,wi)
          v2=0
          do j=1,m_tnk(ii,wi)
            v2=v2+phi(jj,j)*km(j,i)
          end do
          v1=v1+v2*phi(jj,i)
          !print*,v2,tmp_v(jj,i)
        end do
        !vcva(wi,mm+jj,mm+jj)=v1
        vcva(wi,mm+jj)=v1
        !print*,v1,tmp_vcva(jj,jj)
      end do
      !$OMP END PARALLEL DO

      !deallocate(tmp_v) !,tmp_vcva2)
      !deallocate(phi,km,tmp_wk,tmp_vcva)
      deallocate(phi,km,tmp_wk)
    end do

    !This is not used in univariate, for multivariate need to update >>> check
    do ii=2,ttn
      do jj=1,ii-1
        !see aireml_mrnm7.f90 
      end do
    end do     !this is not used at the moment >>> check
  end do  !wi
  call cpu_time (t1_cpu)
  print*,'pre-process done - time:',real(t1_cpu-t2_cpu) !,LA
  !pause

  call cpu_time (t2_cpu)
  !pm=0;vmi=0;trnx=0
  vmi=0;trnx=0
  do xi=1,tn
    trnx=trnx+rnm(xi)
    trny=0
    do yi=1,xi
      trny=trny+rnm(yi)
      !$OMP PARALLEL DO PRIVATE(i) 
      do i=1,rnm(xi)
        !pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)+vcva(1,trnx-rnm(xi)+i,trnx-rnm(xi)+i)+vcve(xi,xi)
        vmi(trnx-rnm(xi)+i)=vmi(trnx-rnm(xi)+i)+vcva(1,trnx-rnm(xi)+i)+vcve(xi,xi)
      end do
      !$OMP END PARALLEL DO
    end do
  end do
  call cpu_time (t1_cpu)
  print*,'V matrix done - time:',real(t1_cpu-t2_cpu) !,LA

  call cpu_time (t2_cpu)
  !call cholesky_inv (pm,trn,LA)
  LA=0
  do i=1,trn
    LA=LA+log(vmi(i))
    !pm(i,i)=1/vmi(i)
    vmi(i)=1/vmi(i)
  end do
  call cpu_time (t1_cpu)
  print*,'V inverse done - time:',real(t1_cpu-t2_cpu) !,LA
  !print*,pm(1,:) !,LA
  !pause
 
    ! P = (VI)-(VI X (X'VI X)I' X' VI)******************************
    !xvm=MATMUL(transpose(xm),pm)
    !$OMP PARALLEL DO PRIVATE (i,j) 
    do i=1,tfn
      do j=1,trn
        !v1=0
        !do k=1,trn
        !  v1=v1+xm(k,i)*pm(k,j)    !xm*pm
        !end do
        xvm(i,j)=xm(j,i)*vmi(j)
      end do
    end do
    !$OMP END PARALLEL DO

    !xvmx=MATMUL(xvm,xm)
    !$OMP PARALLEL DO PRIVATE (i,j,k,v1) 
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

    call cholesky_inv (xvmx,tfn,LG)
    !BLUE = (X'VI X)I XVI y **********************
    !beta=matmul(matmul(matmul(xvmx,transpose(xm)),pm),tyv)
    !$OMP PARALLEL DO PRIVATE (i,j,k,v1) 
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

    !$OMP PARALLEL DO PRIVATE (v2) 
    do i=1,tfn
      v2=0
      do j=1,trn
        !v1=0
        !do k=1,trn
        !  v1=v1+xvm2(i,k)*pm(k,j)    !*pm
        !end do
        !v2=v2+v1*tyv(j,1)           !*tyv
        v2=v2+xvm2(i,j)*vmi(j)*tyv(j,1)           !*tyv
        xvm3(j,i)=xvm2(i,j)*vmi(j)
        !print*,xvm3(j,i)
      end do
      beta(i,1)=v2
    end do
    !$OMP END PARALLEL DO
   !*********************************************
   !print*,'beta:',beta

    v2=0
    res=tyv-matmul(xm,beta)
    tr_vmi=sum(vmi)
    do i=1,trn
      py(i,1)=vmi(i)*res(i,1)
      v2=v2+py(i,1)*tyv(i,1)
      !llik(i)=-0.5*(log(1/vmi(i))+py(i,1)*tyv(i,1)) !REML, ignoring ln(XV-1X)
      llik(i)=-0.5*(log(1/vmi(i))+LG/trn+py(i,1)*tyv(i,1)) !REML, equal contribution for ln(XV-1X)
      llik2(i)=-0.5*(log(1/vmi(i))+LG*(vmi(i)/tr_vmi)+py(i,1)*tyv(i,1)) !REML, weighted for ln(XV-1X)
      !llik2(i)=-0.5*(log(1/vmi(i))+(res(i,1)*res(i,1))*vmi(i))     !maximum likelihood
    end do
    !print*,py(1:10,1)

    LKH=-0.5*(LA+LG+v2)
    !print*,LKH,LA,LG,v2,log(sum(vmi)),sum(log(1/vmi))  !ypyf
    call cpu_time (t2_cpu)
    !print*,'LKH:',t2_cpu-t1_cpu

    if ((LKH>-99999999 .and. LKH<99999999) .and. (LKH>=LKHP)) then
      itit=1
      !print*,'ok',itit,LKHP
    else
      itit=itit+1
      LKH=LKHP
      print*,itit-1,'likelihood NaN >> update reduced by the factor'
      !print*,'not ok',itit,LKHP
      goto 111
    end if

    PRINT '(a7,100f14.4)','LKH',LKH
    do xi=1,tn
      PRINT '(a7,100f14.4)','Ve',vcve(xi,xi)
    end do
    print*,''
    do wi=1,mn
      do ii=1,ttn
        do xi=1,m_tnk(ii,wi)
          PRINT '(a7,100f14.4)','Vk',m_tkm(wi,ii,xi,xi)
        end do
        do xi=2,m_tnk(ii,wi)
          do yi=1,xi-1
            PRINT '(a7,100f14.4)','cov',m_tkm(wi,ii,xi,yi)
          end do
        end do
        print*,''
      end do
      do ii=2,ttn
        do jj=1,ii-1
          do xi=1,m_tnk(ii,wi)
            do yi=1,m_tnk(jj,wi)
              PRINT '(a7,100f14.4)','CVk',ctkm(wi,ii,jj,xi,yi)
            end do
          end do
          print*,''
        end do
      end do
      print*,''
    end do

    ! AI matrix

    !making projection matrix (P) from V
    !$OMP PARALLEL DO PRIVATE (i,j,k,v1) 
    !do i=1,trn
      do j=1,trn
        v1=0
        do k=1,tfn
          v1=v1+xvm3(j,k)*xvm(k,j)
        end do
        !pm(i,j)=pm(i,j)-v1
        pm2(j)=vmi(j)-v1
        !print*,v1
      end do
      !pm2(i)=pm(i,i)
    !end do
    !$OMP END PARALLEL DO

    !py=MATMUL(pm,tyv)       !for py
    !print*,py(1:10,1)

    !print*,'ai matrix'
    call cpu_time (t1_cpu)
    aim=0
    dldv=0

    !Ve (diagonal)  **********************************
    !simplification (from aireml_m_eig4.f90)
    !making dPy (see thompson3.xlsx) 
    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)
      dpy=0
      xvmdpy=0  !xvm %*% dpy
      do i=1,rnm(xi)
        dpy(trnx-rnm(xi)+i,1)=py(trnx-rnm(xi)+i,1)
        v1=0
        do j=1,tfn
          xvmdpy(j,1)=xvmdpy(j,1)+xvm(j,trnx-rnm(xi)+i)*dpy(trnx-rnm(xi)+i,1)
        end do
      end do
      d_beta=matmul(xvmx,xvmdpy)   !beta with dpy as response variable
      d_res=dpy-matmul(xm,d_beta)  !res with dpy as repsonse variable

      !making PdPy
      v3=0
      do ui=1,tn
        do i=1,trn    !same as pedn
          v2=0
          do vi=1,tn
            v2=v2+vmi(i)*d_res(vi*pedn-pedn+i,1)
          end do
          pdpy(ui*pedn-pedn+i,1)=v2
          v3=v3+dpy(ui*pedn-pedn+i,1)*v2
        end do
      end do
      aim(xi,xi)=v3
      !print*,v3
    end do

    !dldv Ve ***********************************************************
    trnx=0
    do xi=1,tn            !tn strictly 1, check >>>
      trnx=trnx+rnm(xi)   !total rn upto xi
      tr1=0;v2=0
      do i=1,rnm(xi)
        !tr1=tr1+pm2(trnx-rnm(xi)+i)
        !v2=v2+py(trnx-rnm(xi)+i,1)*py(trnx-rnm(xi)+i,1)
        dldv_sc(i,xi)=0.5*(py(trnx-rnm(xi)+i,1)*py(trnx-rnm(xi)+i,1)-pm2(trnx-rnm(xi)+i))
      end do 
      !tr1=-0.5*tr1
      dldv(xi,1)=sum(dldv_sc(:,xi))    !v2*(0.5)+tr1
      !print*,v2,tr1,'dldv:',dldv(xi,1),sum(dldv_sc(:,xi))
    end do

    call cpu_time (t2_cpu)
    print*,'derivatives for ve - time:',real(t2_cpu-t1_cpu)

    !for Vg
    nb2=tn
    do k=1,mn !no. random effects
      !print*,'!for Vg (doagonal) ****************************************** '

      do ii=1,ttn
        allocate(phi(rnm(ii),m_tnk(ii,k)),km2(m_tnk(ii,k),m_tnk(ii,k)))
        allocate(tmp_wk(rnm(ii)))

        tmp_wk=m_wk(ii,1:rnm(ii))
        call pol (m_tnk(ii,k),tmp_wk,rnm(ii),phi)

        do xi=1,m_tnk(ii,k)   !check
          nb2=nb2+1

          km2=0;km2(xi,xi)=1;pig_vcva=0

          mm=sum(rnm(1:ii))-rnm(ii)
          do i=1,rnm(ii)
            v1=phi(i,xi)*phi(i,xi)
            pig_vcva(k,mm+i)=v1
          end do

          pig=0;trnx=0
          do ui=1,tn
            trnx=trnx+rnm(ui)
              do i=1,rnm(ui)
                pig(trnx-rnm(ui)+i)=pig(trnx-rnm(ui)+i)+pig_vcva(1,trnx-rnm(ui)+i)
              end do
          end do

          idum=tn
          do ll=1,mn
            do kk=1,ttn
              if (ll==k .and. kk==ii) goto 50
              idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
            end do
            do kk=2,ttn
              do mm=1,kk-1
                idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
              end do
            end do
          end do
50        continue

        dpy=0
        xvmdpy=0    !xvm %*% dpy
        do i=1,trn
          !eig simplification
          v10=py(i,1)*pig_vcva(1,i)     !py*pig
          dpy(i,1)=v10
          !xvm %*% dpy
          do j=1,tfn
            xvmdpy(j,1)=xvmdpy(j,1)+xvm(j,i)*dpy(i,1)
          end do
        end do
        d_beta=matmul(xvmx,xvmdpy)   !beta with dpy as response variable
        d_res=dpy-matmul(xm,d_beta)  !res with dpy as repsonse variable

        !making PdPy
        v3=0
        do ui=1,tn
          do i=1,trn    !same as pedn
            v2=0
            do vi=1,tn
              v2=v2+vmi(i)*d_res(i,1)
            end do
            pdpy(i,1)=v2
            v3=v3+dpy(i,1)*v2
          end do
        end do
        !print*,v3
        aim(idum+xi,idum+xi)=v3


          !print*,'vg',idum+xi,tn+nb*k-nb+xi  !'tmpm:',tmpm(1,1)/2

          !for Vg x Ve (off diagonal) !************************'
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnv-rnm(vi)+i,1) !*py
            end do
            aim(idum+xi,vi)=v2
            !print*,'vg x ve:',idum+xi,vi, v2
          end do !vi

          !for Vg(1~[i-1]) x Vg(i) within trait ***'
          do vi=1,xi-1
            km2=0;km2(vi,vi)=1;pig_vcva=0

            mm=sum(rnm(1:ii))-rnm(ii)
            do i=1,rnm(ii)
              v1=phi(i,vi)*phi(i,vi)
              pig_vcva(k,mm+i)=v1
            end do

            v2=0
            do i=1,trn !rnm(vi)
              !eig simplification  >>> check when using multiple GRMs
              v1=pdpy(i,1)*pig_vcva(k,i)   !*pig_tmp
              v2=v2+v1*py(i,1)                                !*py
            end do
              aim(idum+xi,idum+vi)=v2
              !print*,'vg x vg:',v2
            end do

          !dldv Vg *********************************************************
          tr1=0
          do i=1,trn
            !tr1=tr1+pm2(i)*pig(i)
            dldv_sc(i,idum+xi)=0.5*(py(i,1)*pig(i)*py(i,1)-pm2(i)*pig(i))
          end do
          !tr1=-0.5*tr1

          !v1=0
          !do i=1,trn
          !  v1=v1+py(i,1)*pig(i)*py(i,1)
          !end do
        
          dldv(idum+xi,1)=sum(dldv_sc(:,idum+xi))   !v1*0.5+tr1 !************************************
          !print*,'dldv:',dldv(idum+xi,1),sum(dldv_sc(:,idum+xi))
        end do !xi


        !for cov (diagonal) *************************************************
        nbc=0  !for no. block for covariance, i.e. 1: (2,1), 2: (3,1), 3:(3,2),... 
        do xi=2,m_tnk(ii,k)
          do yi=1,xi-1
            nbc=nbc+1
            km2=0;km2(xi,yi)=1;km2(yi,xi)=1;pig_vcva=0

            mm=sum(rnm(1:ii))-rnm(ii)
            do i=1,rnm(ii)
              v1=2*phi(i,xi)*phi(i,yi)
              pig_vcva(k,mm+i)=v1
            end do

          pig=0;trnx=0
          do ui=1,tn
            trnx=trnx+rnm(ui)
            do i=1,rnm(ui)
              pig(trnx-rnm(ui)+i)=pig(trnx-rnm(ui)+i)+pig_vcva(1,trnx-rnm(ui)+i)
            end do
          end do

            idum=tn
            do ll=1,mn
              do kk=1,ttn
                if (ll==k .and. kk==ii) goto 53
                idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
              end do
              do kk=2,ttn
                do mm=1,kk-1
                  idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
                end do
              end do
            end do
53          continue

            dpy=0
            xvmdpy=0    !xvm %*% dpy
            do i=1,trn
              !eig simplification
              v10=py(i,1)*pig_vcva(1,i)     !py*pig
              dpy(i,1)=v10
              !xvm %*% dpy
              do j=1,tfn
                xvmdpy(j,1)=xvmdpy(j,1)+xvm(j,i)*dpy(i,1)
              end do
            end do
            d_beta=matmul(xvmx,xvmdpy)   !beta with dpy as response variable
            d_res=dpy-matmul(xm,d_beta)  !res with dpy as repsonse variable

            !making PdPy
            v3=0
            do ui=1,tn
              do i=1,trn    !same as pedn
                v2=0
                do vi=1,tn
                  v2=v2+vmi(i)*d_res(i,1)
                end do
                pdpy(i,1)=v2
                v3=v3+dpy(i,1)*v2
              end do
            end do
            aim(idum+m_tnk(ii,k)+nbc,idum+m_tnk(ii,k)+nbc)=v3     !check
            !print*,'cov***:',v3

            !for cov x Ve (off diagonal) ************************************
            trnv=0
            do vi=1,tn
              trnv=trnv+rnm(vi)

              v2=0
              do i=1,rnm(vi)
                v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnv-rnm(vi)+i,1) !*py
              end do
              aim(idum+m_tnk(ii,k)+nbc,vi)=v2          !check with Ve1
              !print*,'cov x ve:',v2
            end do

            !for cov x Vg **************************************************
            do vi=1,m_tnk(ii,k)

              km2=0;km2(vi,vi)=1;pig_vcva=0

              mm=sum(rnm(1:ii))-rnm(ii)
              do i=1,rnm(ii)
                v1=phi(i,vi)*phi(i,vi)
                pig_vcva(k,mm+i)=v1
              end do


            v2=0
            do i=1,trn
              !eig simplification  >>> check when using multiple GRMs
              v1=pdpy(i,1)*pig_vcva(k,i)   !*pig_tmp
              v2=v2+v1*py(i,1)                                !*py
            end do

              aim(idum+m_tnk(ii,k)+nbc,idum+vi)=v2  !with the same Vg1 for first trait
              !print*,'cov x vg:',v2,vi,k,pig_vcva(1,1:10)
            end do

            !for cov x cov *************************************************
            !trnv=rnm(1)
            nbc2=0    !no block for cov (vi,ui)
            do vi=2,m_tnk(ii,k)
              !trnv=trnv+rnm(vi)
              trnu=0
              do ui=1,(vi-1)
                nbc2=nbc2+1

                if (nbc2<nbc) then  !off diagonal within cova ************************
                  trnu=trnu+rnm(ui)

                  km2=0;km2(vi,ui)=1;km2(ui,vi)=1;pig_vcva=0

                  mm=sum(rnm(1:ii))-rnm(ii)
                  do i=1,rnm(ii)
                    v1=2*phi(i,vi)*phi(i,ui)
                    pig_vcva(k,mm+i)=v1
                  end do

                  pig_tmp=0;trnx=0
                  do wj=1,tn
                    trnx=trnx+rnm(wj)
                    do i=1,rnm(wj)
                      pig_tmp(trnx-rnm(wj)+i)=pig_tmp(trnx-rnm(wj)+i)+pig_vcva(1,trnx-rnm(wj)+i)
                    end do
                  end do

                  v1=0
                  do i=1,trn
                    !v1=v1+py(i,1)*pig(i,i)*pm(i,i)*pig_tmp(i,i)*py(i,1)
                    v1=v1+py(i,1)*pig(i)*pm2(i)*pig_tmp(i)*py(i,1)
                  end do
                  tmpm(1,1)=v1

                  aim(idum+m_tnk(ii,k)+nbc,idum+m_tnk(ii,k)+nbc2)=tmpm(1,1) !with the wi th cov for first trait
                  !print*,'cov x cov:'
                end if ! *************************************************************
              end do !ui
            end do !vi

          ! dldv for cov
          !tr1=0
          do i=1,trn
            !tr1=tr1+pm2(i)*pig(i)
            dldv_sc(i,idum+m_tnk(ii,k)+nbc)=0.5*(py(i,1)*pig(i)*py(i,1)-pm2(i)*pig(i))
          end do
          !tr1=-0.5*tr1

          !v1=0
          !do i=1,trn
          !  v1=v1+py(i,1)*pig(i)*py(i,1)
          !end do
          dldv(idum+m_tnk(ii,k)+nbc,1)=sum(dldv_sc(:,idum+m_tnk(ii,k)+nbc))   !v1*0.5+tr1
          !print*,dldv(idum+m_tnk(ii,k)+nbc,1),sum(dldv_sc(:,idum+m_tnk(ii,k)+nbc))
        end do !yi
      end do !xi

      deallocate(phi,km2,tmp_wk)

    end do !ii

    end do !k

    call cpu_time (t2_cpu)
    print*,'derivatives done - time:',real(t2_cpu-t1_cpu)
    print*,''

    idum=tn
    do ll=1,mn
      do kk=1,ttn
        idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
      end do
      do kk=2,ttn
        do mm=1,kk-1
          !idum=idum+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
          idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
        end do
      end do
    end do
 
    do i=1,idum
      aim(i,i)=aim(i,i)/2
      do j=1,i-1
        aim(i,j)=aim(i,j)/2
        aim(j,i)=aim(i,j)
      end do
    end do

    !do i=1,idum
    !  print*,i,real(aim(i,1:i))
    !end do
    !print*,''

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
      do ii=1,ttn
        do xi=1,m_tnk(ii,wi)
          wj=wj+1
          if (fix_m_tkm(wi,ii,xi,xi).ne.-99999) then
            aim(:,wj)=0
            aim(wj,:)=0
            aim(wj,wj)=1
          end if
        end do
        do xi=2,m_tnk(ii,wi)
          do yi=1,xi-1
            wj=wj+1
            if (fix_m_tkm(wi,ii,xi,yi).ne.-99999) then
              aim(:,wj)=0
              aim(wj,:)=0
              aim(wj,wj)=1
            end if
          end do
        end do
      end do !ii
      do ii=2,ttn
        do jj=1,ii-1
          do xi=1,m_tnk(ii,wi)
            do xi2=1,m_tnk(jj,wi)
              wj=wj+1
              if (fix_ctkm(wi,ii,jj,xi,xi2).ne.-99999) then
                aim(:,wj)=0
                aim(wj,:)=0
                aim(wj,wj)=1
              end if
            end do
          end do
        end do
      end do 
    end do


    !do i=1,idum
    !  print*,i,real(aim(i,1:i))
    !end do
    !print*,''

    infm=aim
    call cholesky_inv (aim,idum,x1)
    print *,''

    if (LKH-MLKH.gt.0) then
      MLKH=LKH
      sdmI=aim
    end if

    !if (ABS(LKH-LKHP).lt.0.001) goto 1000
    if (ABS(LKH-LKHP).lt.conv) goto 1000


    up=MATMUL(aim,dldv)
    !do i=1,idum
    !  print*,real(up(i,1)),real(dldv(i,1))
    !end do
    !print*,''
    !pause

111 continue
    if (itit==1) then
      vcve_rsv=vcve
      !tkm_rsv=tkm
      m_tkm_rsv=m_tkm
      ctkm_rsv=ctkm
      if (zi.ge.nit) then
        print*,'Likelihood not converged >>> may need a longer iterations'
        print*,''
        goto 1000
      end if

    else
      up=up*0.7**(itit-1)
      vcve=vcve_rsv
      !tkm=tkm_rsv
      m_tkm=m_tkm_rsv
      ctkm=ctkm_rsv
    end if

    do xi=1,tn
      vcve(xi,xi)=vcve(xi,xi)+up(xi,1)
      !if (vcve(xi,xi)<0) vcve(xi,xi)=0.000001
    end do
    do wi=1,mn
      do ii=1,ttn

        idum=tn
        do ll=1,mn
          do kk=1,ttn
            if (ll==wi .and. kk==ii) goto 80
            idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
          end do
          do kk=2,ttn
            do mm=1,kk-1
              idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
            end do
          end do
        end do
80      continue

        do xi=1,m_tnk(ii,wi)
          m_tkm(wi,ii,xi,xi)=m_tkm(wi,ii,xi,xi)+up(idum+xi,1)
          !print*,idum+xi
        end do
        nbc=0
        do xi=2,m_tnk(ii,wi)
          do yi=1,xi-1
            nbc=nbc+1
            m_tkm(wi,ii,xi,yi)=m_tkm(wi,ii,xi,yi)+up(idum+m_tnk(ii,wi)+nbc,1)
            m_tkm(wi,ii,yi,xi)=m_tkm(wi,ii,xi,yi)
            !print*,idum+m_tnk(ii,wi)+nbc
          end do
        end do
      end do ! ii

      do ii=1,ttn
        do jj=1,ii-1
          nbc2=0
          do xi=1,m_tnk(ii,wi)
            do yi=1,m_tnk(jj,wi)
              nbc2=nbc2+1

              idum2=tn
              do ll=1,mn
                do kk=1,ttn
                   idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                end do
                do kk=2,ttn
                  do mm=1,kk-1
                    if (ll==wi .and. kk==ii .and. mm==jj) goto 81
                    idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                  end do
                end do
              end do
81            continue

              ctkm(wi,ii,jj,xi,yi)=ctkm(wi,ii,jj,xi,yi)+up(idum2+nbc2,1)
              !print*,'*',idum+m_tnk(ii,wi)+nbc+nbc2,idum2+nbc2
            end do
          end do

        end do ! jj
      end do ! ii
    end do !wi
  END do !zi


  !PRINT*,'LKH not converged na na na'
  !WRITE(41,*)'LKH not converged na na na'
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
    do ii=1,ttn

          idum=tn
          do ll=1,mn
            do kk=1,ttn
              if (ll==wi .and. kk==ii) goto 95
              idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
            end do
            do kk=2,ttn
              do mm=1,kk-1
                idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
              end do
            end do
          end do
95        continue

      do xi=1,m_tnk(ii,wi)
        write(41,'(a7,100f14.4)')'Vk',m_tkm(wi,ii,xi,xi),sqrt(sdmI(idum+xi,idum+xi))
      end do

      nbc=0
      do xi=2,m_tnk(ii,wi)
        do yi=1,xi-1
          nbc=nbc+1
          write(41,'(a7,100f14.4)')'cov',m_tkm(wi,ii,xi,yi),sqrt(sdmI(idum+m_tnk(ii,wi)+nbc,idum+m_tnk(ii,wi)+nbc))
        end do
      end do

      write(41,*)''
    end do ! ii

    do ii=2,ttn
      do jj=1,ii-1
        nbc2=0
        do xi=1,m_tnk(ii,wi)
          do yi=1,m_tnk(jj,wi)

            idum2=tn
            do ll=1,mn
              do kk=1,ttn
                 idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
              end do

              do kk=2,ttn
                do mm=1,kk-1
                  if (ll==wi .and. kk==ii .and. mm==jj) goto 96
                  idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                end do
              end do
            end do
96          continue

            nbc2=nbc2+1
            write(41,'(a7,100f14.4)')'CVk',ctkm(wi,ii,jj,xi,yi),sqrt(sdmI(idum2+nbc2,idum2+nbc2))
          end do
        end do
      end do !jj
      write(41,*)''
    end do !ii

  end do ! wi

  write(41,'(a7,100f14.4)')'LKH',LKH !,real(zi) !indicate convergence
  write(41,*)''



  idum=tn
  do ll=1,mn
    do kk=1,ttn
      idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
    end do
    do kk=2,ttn
      do mm=1,kk-1
        idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
      end do
    end do
  end do

  write(41,*)''
  do i=1,idum
    write(41,'(100g18.8)')(sdmI(i,1:i))
  end do

  close(41)


  open (unit=55,file=trim(fl5)//".vuong_2.1",status='unknown')
    write(55,'(100g18.8)')infm(1,1),infm(1,3:4)
    write(55,'(100g18.8)')infm(3,1),infm(3,3:4)
    write(55,'(100g18.8)')infm(4,1),infm(4,3:4)
  close(55)

  open (unit=56,file=trim(fl5)//".vuong_2.2",status='unknown')
  do i=1,rn
    write(56,'(100g18.8)')dldv_sc(i,1),dldv_sc(i,3:4)
  end do
  close(56)

  open (unit=57,file=trim(fl5)//".vuong_llik",status='unknown')
  do i=1,rn
    write(57,'(100g18.8)')llik(i)
  end do
  close(57)

  open (unit=58,file=trim(fl5)//".vuong_llik2",status='unknown')
  do i=1,rn
    write(58,'(100g18.8)')llik2(i)
  end do
  close(58)


if (fl11.ne.'null') then
  print*,'for BLUP solution after convergence *******************'
  !!make it triangle -> full ****************
  !do xi=1,tn
  !  do yi=1,xi
  !    do k=1,mn
  !      vcva(k,yi,xi)=vcva(k,xi,yi)
  !      !print*,vcva
  !    end do
  !  end do
  !end do !***********************************

  do k=1,mn

    !allocate(gr_km_phi(sum(m_tnk(1:ttn,k)),tn))
    allocate(gr_km_phi(sum(m_tnk(1:ttn,k)),trn))

    do ii=1,ttn
      allocate(phi(rnm(ii),m_tnk(ii,k)),km(m_tnk(ii,k),m_tnk(ii,k)))
      allocate(tmp_wk(rnm(ii)),km_phi(m_tnk(ii,k),rnm(ii)))

      tmp_wk=m_wk(ii,1:rnm(ii))
      call pol (m_tnk(ii,k),tmp_wk,rnm(ii),phi)
      km=m_tkm(k,ii,1:m_tnk(ii,k),1:m_tnk(ii,k))

      km_phi=matmul(km,transpose(phi))

      do i=1,m_tnk(ii,k)
        do j=1,rnm(ii)
          i2=sum(m_tnk(1:ii,k))-m_tnk(ii,k)+i
          j2=sum(rnm(1:ii))-rnm(ii)+j
          gr_km_phi(i2,j2)=km_phi(i,j)
        end do
      end do
      deallocate(phi,km,tmp_wk,km_phi)
    end do ! ii

    do ii=2,ttn   !two traits only, check when > 2 trait 
      do jj=1,ii-1
        allocate(phi(rnm(jj),m_tnk(jj,k)),km(m_tnk(ii,k),m_tnk(jj,k)))
        allocate(tmp_wk(rnm(jj)),km_phi(m_tnk(ii,k),rnm(jj)))

        tmp_wk=m_wk(jj,1:rnm(jj))
        call pol (m_tnk(jj,k),tmp_wk,rnm(jj),phi)

        allocate(phi2(rnm(ii),m_tnk(ii,k)),tmp_wk2(rnm(ii)))
        allocate(km_phi2(m_tnk(jj,k),rnm(ii)))

        tmp_wk2=m_wk(ii,1:rnm(ii))
        call pol (m_tnk(ii,k),tmp_wk2,rnm(ii),phi2)

        km=ctkm(k,ii,jj,1:m_tnk(ii,k),1:m_tnk(jj,k))

        km_phi=matmul(km,transpose(phi))

        do i=1,m_tnk(ii,k)
          do j=1,rnm(jj)
            i2=sum(m_tnk(1:ii,k))-m_tnk(ii,k)+i
            j2=sum(rnm(1:jj))-rnm(jj)+j

            gr_km_phi(i2,j2)=km_phi(i,j)
          end do
        end do

        km_phi2=matmul(transpose(km),transpose(phi2))

        do i=1,m_tnk(jj,k)
          do j=1,rnm(ii)
            i2=sum(m_tnk(1:jj,k))-m_tnk(jj,k)+i
            j2=sum(rnm(1:ii))-rnm(ii)+j
            gr_km_phi(i2,j2)=km_phi2(i,j)
          end do
        end do

        deallocate(phi,km,tmp_wk,km_phi,phi2,tmp_wk2,km_phi2)

      end do ! ii
    end do ! jj

    do xi=1,sum(m_tnk(1:ttn,k))    !check
      do zi=1,pedn
        v1=0;v2=0

        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          !do j=1,rnm(yi)
          !  !v1=v1+mbin(k,yidx(yi,j),zi)*gr_km_phi(xi,yi)*py(trny-rnm(yi)+j,1)
          !  !v1=v1+mbin(k,yidx(yi,j),zi)*gr_km_phi(xi,j)*py(trny-rnm(yi)+j,1)
          !  v1=v1+mbin(k,yidx(yi,j),zi)*gr_km_phi(xi,trny-rnm(yi)+j)*py(trny-rnm(yi)+j,1)
          !end do
          v1=v1+gr_km_phi(xi,trny-rnm(yi)+zi)*py(trny-rnm(yi)+zi,1)
        end do
        blup_ebv(zi,xi,k)=v1
      end do !zi
    end do !xi
    deallocate(gr_km_phi)

  end do !k

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
  write(43,'(a6,a14,a9,a14,a24)')'trait','random eff.#','RR order','ordered ind#','EBVs'

  do ii=1,ttn
    do k=1,mn
      do xi=1,m_tnk(ii,k)    !check
        do i=1,pedn
          write(43,'(i6,i14,i9,i14,f24.16)')ii,k,xi,i,blup_ebv(i,sum(m_tnk(1:ii,k))-m_tnk(ii,k)+xi,k)
        end do
      end do
    end do
  end do
  close(43)

end if

if (fl11r.ne.'null') then
  print*,'for BLUP solution and reliability after convergence *************'
  print*,'this should be obtained from MRNM *************'

  do k=1,mn
    allocate(gr_km_phi(sum(m_tnk(1:ttn,k)),trn))
    do ii=1,ttn
      allocate(phi(rnm(ii),m_tnk(ii,k)),km(m_tnk(ii,k),m_tnk(ii,k)))
      allocate(tmp_wk(rnm(ii)),km_phi(m_tnk(ii,k),rnm(ii)))

      tmp_wk=m_wk(ii,1:rnm(ii))
      call pol (m_tnk(ii,k),tmp_wk,rnm(ii),phi)
      km=m_tkm(k,ii,1:m_tnk(ii,k),1:m_tnk(ii,k))

      km_phi=matmul(km,transpose(phi))

      do i=1,m_tnk(ii,k)
        do j=1,rnm(ii)
          i2=sum(m_tnk(1:ii,k))-m_tnk(ii,k)+i
          j2=sum(rnm(1:ii))-rnm(ii)+j
          gr_km_phi(i2,j2)=km_phi(i,j)
        end do
      end do
      deallocate(phi,km,tmp_wk,km_phi)
    end do ! ii

    do ii=2,ttn   !two traits only, check when > 2 trait 
      do jj=1,ii-1
        allocate(phi(rnm(jj),m_tnk(jj,k)),km(m_tnk(ii,k),m_tnk(jj,k)))
        allocate(tmp_wk(rnm(jj)),km_phi(m_tnk(ii,k),rnm(jj)))

        tmp_wk=m_wk(jj,1:rnm(jj))
        call pol (m_tnk(jj,k),tmp_wk,rnm(jj),phi)

        allocate(phi2(rnm(ii),m_tnk(ii,k)),tmp_wk2(rnm(ii)))
        allocate(km_phi2(m_tnk(jj,k),rnm(ii)))

        tmp_wk2=m_wk(ii,1:rnm(ii))
        call pol (m_tnk(ii,k),tmp_wk2,rnm(ii),phi2)

        km=ctkm(k,ii,jj,1:m_tnk(ii,k),1:m_tnk(jj,k))

        km_phi=matmul(km,transpose(phi))

        do i=1,m_tnk(ii,k)
          do j=1,rnm(jj)
            i2=sum(m_tnk(1:ii,k))-m_tnk(ii,k)+i
            j2=sum(rnm(1:jj))-rnm(jj)+j

            gr_km_phi(i2,j2)=km_phi(i,j)
          end do
        end do

        km_phi2=matmul(transpose(km),transpose(phi2))

        do i=1,m_tnk(jj,k)
          do j=1,rnm(ii)
            i2=sum(m_tnk(1:jj,k))-m_tnk(jj,k)+i
            j2=sum(rnm(1:ii))-rnm(ii)+j
            gr_km_phi(i2,j2)=km_phi2(i,j)
          end do
        end do

        deallocate(phi,km,tmp_wk,km_phi,phi2,tmp_wk2,km_phi2)

      end do ! ii
    end do ! jj

    do xi=1,sum(m_tnk(1:ttn,k))    !check
      do zi=1,pedn
        v1=0;v2=0

        trny=0
        do yi=1,tn
          !trny=trny+rnm(yi)
          !do j=1,rnm(yi)
          !  v1=v1+mbin(k,yidx(yi,j),zi)*gr_km_phi(xi,trny-rnm(yi)+j)*py(trny-rnm(yi)+j,1)
          !end do
          v1=v1+gr_km_phi(xi,trny-rnm(yi)+zi)*py(trny-rnm(yi)+zi,1)
        end do
        blup_ebv(zi,xi,k)=v1
      end do !zi
    end do !xi

    !to get reliability for BLUP, i.e. GPA
    do xi=1,sum(m_tnk(1:ttn,k))    !check
      do zi=1,pedn
        v2=0
        trnv=0
        do vi=1,tn
          trnv=trnv+rnm(vi)
          !do i=1,rnm(vi)
            v1=0
            trny=0
            do yi=1,tn
              trny=trny+rnm(yi)
              !do j=1,rnm(yi)
              !  !v1=v1+mbin(k,yidx(yi,j),zi)*gr_km_phi(xi,trny-rnm(yi)+j)*py(trny-rnm(yi)+j,1)
              !  v1=v1+mbin(k,yidx(yi,j),zi)*gr_km_phi(xi,trny-rnm(yi)+j)*pm(trny-rnm(yi)+j,trnv-rnm(vi)+i)
              !end do
              !v1=v1+gr_km_phi(xi,trny-rnm(yi)+zi)*pm(trny-rnm(yi)+zi,trnv-rnm(vi)+zi)
              v1=v1+gr_km_phi(xi,trny-rnm(yi)+zi)*pm2(trny-rnm(yi)+zi)   !check if tn>1
            end do
            !v2=v2+v1*mbin(k,yidx(vi,i),zi)
          !end do
          v2=v2+v1
        end do
        blup_r(zi,xi,k)=v2
      end do !zi
    end do !xi

    deallocate(gr_km_phi)

  end do !k

  open (unit=47,file=trim(fl11r)//".fsl",status='unknown')
  do xi=1,tfn
    write(47,'(i3,100f24.16)')xi,beta(xi,1),sqrt(xvmx(xi,xi))
  end do
  close(47)

  open (unit=43,file=fl11r,status='unknown')
  write(43,'(100f24.16)')beta

  do ii=1,ttn
    do k=1,mn
      do xi=1,m_tnk(ii,k)    !check
        do i=1,pedn
          write(43,'(i3,i3,i3,i7,100f24.16)')ii,k,xi,i,blup_ebv(i,sum(m_tnk(1:ii,k))-m_tnk(ii,k)+xi,k)
        end do
      end do
    end do
  end do
  close(43)

  open (unit=46,file=trim(fl11r)//".r2",status='unknown')
  do ii=1,ttn
    do k=1,mn
      do xi=1,m_tnk(ii,k)    !check
        do i=1,pedn
          write(46,'(i3,i3,i3,i7,100f24.16)')ii,k,xi,i,blup_r(i,sum(m_tnk(1:ii,k))-m_tnk(ii,k)+xi,k)**2
        end do
      end do
    end do
  end do
  close(46)


end if







end subroutine

