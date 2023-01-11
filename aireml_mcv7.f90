
subroutine aireml_mcv_m_h (mbin,yidx,yv,rn,rnm,trn,pedn,ttn,tn,fnm,tfn,mn,fl4,fl5,fl11,fl12,fl24,xmm,nit,conv,thn)        

!***********************************************************************
!aireml_rnm: reaction norm model 
!S. Hong Lee (2017)

!mbin : matrices bin
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
integer::rn1,rn2,yidx1(rn),yidx2(rn),ll,mm,mm2,thn
real::mbin(mn,pedn,pedn),yv(rn,tn),tyv(trn,1)

double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),conv
double precision::vcva(mn,trn,trn),pig_vcva(mn,trn,trn),vcve(tn,tn),fix_vcve(tn,tn)
double precision::vcve_rsv(tn,tn)
double precision::blup_ebv(pedn,100,mn),beta(tfn,1) ! check assuming that total sum k (pol oreder) < 100
integer::i,j,m,k,k2,l,io
double precision::sum_h2,sum_v(tn)

double precision::LC,LA,LG,ypy,x1,x2,x3,x10,x11,v1,v2,v10,v11
double precision::y1,y2,y3,y10,y11,z1,z2,z3,z10,z11
double precision::LKHP,MLKH,mva(mn),mva2(mn),mve,mve2,mcov(mn)
double precision::h2(mn),h2n(mn),tr1,tr2(mn),tr3

double precision::xm(trn,tfn),xmm(tn,pedn,10000)

double precision::xvmx(tfn,tfn),xvm(tfn,trn) 
double precision::xvm2(tfn,trn),xvm3(trn,tfn) 
double precision::pm(trn,trn),py(trn,1),py_tmp(1,trn),xb(trn,1)
double precision::py_tmp1(1,trn)
double precision::py_tmp2(1,trn)
double precision::py_tmp3(1,trn)
double precision::pig(trn,trn),pig_tmp(trn,trn)

double precision::v_tmp(trn),v2_tmp(trn)

!time measure **********************************************************
INTEGER::now(8)
CHARACTER*8 date
CHARACTER*20 time,zone
double precision:: timef, t1_real, t2_real, elapsed_real_secs
double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
!***********************************************************************

character(len=7)::cdum
character(len=64)::fl4,fl5,fl11,fl12,fl24
logical::file_exist

double precision::m_tkm(mn,ttn,20,20),m_tkm_rsv(mn,ttn,20,20),fix_m_tkm(mn,ttn,20,20)
double precision::ctkm(mn,ttn,ttn,20,20),ctkm_rsv(mn,ttn,ttn,20,20),fix_ctkm(mn,ttn,ttn,20,20)
double precision,allocatable::km(:,:),km2(:,:),phi(:,:),km2_2(:,:),phi_2(:,:)
double precision,allocatable::phi2(:,:),tmp_vcva(:,:),tmp_vcva_2(:,:),phi_1(:,:),tmp_v(:,:)
!double precision::wk(tn),xk(tn)
double precision::wk(pedn),m_wk(mn,ttn,pedn)
double precision,allocatable::tmp_wk(:),tmp_wk2(:),tmp_wk_2(:),tmp_wk_1(:)
double precision,allocatable::km_phi(:,:),km_phi2(:,:),gr_km_phi(:,:)

double precision,allocatable::aim(:,:)
double precision,allocatable::sdmI(:,:)
double precision,allocatable::up(:,:),dldv(:,:)

integer::idum,idum2,ii,jj,kk,oo,ii2,jj2,xi2,i2,j2

!CDF
real ( kind = 8 ) bound,mean,pv,q,sd,x
integer ( kind = 4 ) status,which


open (UNIT=46,FILE=fl12,STATUS='old')
read(46,*)(sub_tn(i),i=1,ttn)
do i=1,ttn
  if (sub_tn(i).ne.rnm(i)) then
    print*,'# E variable and # records not matched for the trait',i
    pause
  end if
end do

do i=1,ttn
  read(46,*)m_tnk(i,:)
  k2=0
  !do j=1,sub_tn(i)
  do j=1,rn
    if (yv(j,i).ne.-99999) then
      k2=k2+1
      read(46,*,iostat=io)(m_wk(k,i,k2),k=1,mn)
      if (io.ne.0) then
        print*,'E variable and phenotype not matched >>> check'
        pause
      end if
    elseif (yv(j,i).eq.-99999) then
      read(46,*)cdum
    end if
  end do
  if (k2.ne.sub_tn(i)) then
    print*,'E variable and phenotype not matched',i,k2,sub_tn(i)
    pause
  end if
end do
close(46)

!do i=1,ttn
!  m_wk(i,:)=wk   !check >>> for unbalanced E variable, modify this
!  !print*,m_wk(i,:)
!end do
!pause

!print*,tnk
!print*,wk

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

allocate(aim(idum,idum),sdmI(idum,idum),up(idum,1),dldv(idum,1))

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
  vcve=0;vcva=0;m_tkm=0;ctkm=0   !;km=0
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
      allocate(tmp_wk(rnm(ii)),tmp_vcva(rnm(ii),rnm(ii)))

      tmp_wk=m_wk(wi,ii,1:rnm(ii))
      call pol (m_tnk(ii,wi),tmp_wk,rnm(ii),phi)

      km=m_tkm(wi,ii,1:m_tnk(ii,wi),1:m_tnk(ii,wi))

      !tmp_vcva=matmul(matmul(phi,km),transpose(phi))
      allocate(tmp_v(rnm(ii),m_tnk(ii,wi))) !,tmp_vcva2(rnm(ii),rnm(ii)))
      call dgemm ('N','N',rnm(ii),m_tnk(ii,wi),m_tnk(ii,wi),1.0D0,phi,rnm(ii),km,m_tnk(ii,wi),0.0D0,tmp_v,rnm(ii))
      call dgemm ('N','T',rnm(ii),rnm(ii),m_tnk(ii,wi),1.0D0,tmp_v,rnm(ii),phi,rnm(ii),0.0D0,tmp_vcva,rnm(ii))
      deallocate(tmp_v) !,tmp_vcva2)

      mm=sum(rnm(1:ii))-rnm(ii)
      !$OMP PARALLEL DO PRIVATE(jj, kk)
      do jj=1,rnm(ii)
        do kk=1,rnm(ii)
          vcva(wi,mm+jj,mm+kk)=tmp_vcva(jj,kk)
        end do
      end do
      !$OMP END PARALLEL DO
      deallocate(phi,km,tmp_wk,tmp_vcva)
    end do

    do ii=2,ttn
      do jj=1,ii-1
        allocate(km(m_tnk(ii,wi),m_tnk(jj,wi)),tmp_vcva(rnm(ii),rnm(jj)))
        allocate(phi(rnm(ii),m_tnk(ii,wi)),tmp_wk(rnm(ii)))
        tmp_wk=m_wk(wi,ii,1:rnm(ii))
        call pol (m_tnk(ii,wi),tmp_wk,rnm(ii),phi)

        allocate(phi2(rnm(jj),m_tnk(jj,wi)),tmp_wk2(rnm(jj)))
        tmp_wk2=m_wk(wi,jj,1:rnm(jj))
        call pol (m_tnk(jj,wi),tmp_wk2,rnm(jj),phi2)

        km=ctkm(wi,ii,jj,1:m_tnk(ii,wi),1:m_tnk(jj,wi))

        !tmp_vcva=matmul(matmul(phi,km),transpose(phi2))
        allocate(tmp_v(rnm(ii),m_tnk(jj,wi)))
        call dgemm ('N','N',rnm(ii),m_tnk(jj,wi),m_tnk(ii,wi),1.0D0,phi,rnm(ii),km,m_tnk(ii,wi),0.0D0,tmp_v,rnm(ii))
        call dgemm ('N','T',rnm(ii),rnm(jj),m_tnk(jj,wi),1.0D0,tmp_v,rnm(ii),phi2,rnm(jj),0.0D0,tmp_vcva,rnm(ii))
        deallocate(tmp_v)

        mm=sum(rnm(1:ii))-rnm(ii)
        mm2=sum(rnm(1:jj))-rnm(jj)
        !$OMP PARALLEL DO PRIVATE(ll, kk)
        do ll=1,rnm(ii)
          do kk=1,rnm(jj)
            vcva(wi,mm+ll,mm2+kk)=tmp_vcva(ll,kk)
            vcva(wi,mm2+kk,mm+ll)=tmp_vcva(ll,kk)
          end do
        end do
        !$OMP END PARALLEL DO
        deallocate(phi,phi2,km,tmp_wk,tmp_wk2,tmp_vcva)
      end do
    end do
  end do  !wi
  call cpu_time (t1_cpu)
  print*,'pre-process done - time:',real(t1_cpu-t2_cpu) !,LA
  !pause

  !print*,387
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
            pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)+mbin(k,yidx(xi,i),yidx(yi,j))*vcva(k,trnx-rnm(xi)+i,trny-rnm(yi)+j)
          end do
          if (xi.ne.yi) then
            pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)
          end if
        end do
        if (xi==yi) then
          pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)+vcve(xi,xi)
        end if
      end do
      !$OMP END PARALLEL DO
    end do
  end do
  call cpu_time (t1_cpu)
  print*,'V matrix done - time:',real(t1_cpu-t2_cpu) !,LA
  !print*,418

  call cpu_time (t2_cpu)
  call cholesky_inv (pm,trn,LA)
  call cpu_time (t1_cpu)
  print*,'V inverse done - time:',real(t1_cpu-t2_cpu) !,LA
  !print*,LA
  !pause
 
    ! P = (VI)-(VI X (X'VI X)I' X' VI)******************************
    !xvm=MATMUL(transpose(xm),pm)
    !$OMP PARALLEL DO PRIVATE (i,j,k,v1) 
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

    !$OMP PARALLEL DO PRIVATE (i,j,k,v1) 
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
    !v2=0
    !do i=1,trn
    !  v1=0
    !  do j=1,trn
    !    v1=v1+tyv(j,1)*pm(j,i)
    !  end do
    !  v2=v2+v1*tyv(i,1)
    !end do

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
    !print*,LKH,LA,LG,v2  !ypyf
    !pause
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
    py=MATMUL(pm,tyv)       !for py

    !print*,'ai matrix'
    call cpu_time (t1_cpu)
    aim=0
    dldv=0

    !Ve (diagonal)  **********************************
    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)   !total rn upto xi
      v2=0
      do i=1,rnm(xi)
        v1=0
        do j=1,rnm(xi)
          v1=v1+py(trnx-rnm(xi)+j,1)*pm(trnx-rnm(xi)+j,trnx-rnm(xi)+i)
        end do
        v2=v2+v1*py(trnx-rnm(xi)+i,1)
      end do 
      aim(xi,xi)=v2
      !print*,aim(xi,xi)/2
    end do

    !dldv Ve ***********************************************************
    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)   !total rn upto xi
      tr1=0
      do i=1,rnm(xi)
        tr1=tr1+pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)
      end do
      tr1=-0.5*tr1
      v2=0
      do i=1,rnm(xi)
        v2=v2+py(trnx-rnm(xi)+i,1)*py(trnx-rnm(xi)+i,1)
      end do 
      dldv(xi,1)=v2*(0.5)+tr1
      !print*,v2,tr1 !'dldv:',dldv(xi,1)
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
            v1=v1+py(trny-rnm(yi)+j,1)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
          end do
          v2=v2+v1*py(trnx-rnm(xi)+i,1)
        end do 
        aim(xi,yi)=v2
        !print*,'ve x ve',xi,yi,v2/2
      end do
    end do
    !pause
    call cpu_time (t2_cpu)
    print*,'derivatives for ve - time:',real(t2_cpu-t1_cpu)

    !for Vg
    nb2=tn
    do k=1,mn !no. random effects
      !for Vg (doagonal) ******************************************************* 

      do ii=1,ttn
        allocate(phi(rnm(ii),m_tnk(ii,k)),km2(m_tnk(ii,k),m_tnk(ii,k)))
        allocate(tmp_wk(rnm(ii)),tmp_vcva(rnm(ii),rnm(ii)))

        tmp_wk=m_wk(k,ii,1:rnm(ii))
        call pol (m_tnk(ii,k),tmp_wk,rnm(ii),phi)

        do xi=1,m_tnk(ii,k)   !check
          nb2=nb2+1

          km2=0;km2(xi,xi)=1;pig_vcva=0
          !tmp_vcva=matmul(matmul(phi,km2),transpose(phi))
          allocate(tmp_v(rnm(ii),m_tnk(ii,k)))
          call dgemm ('N','N',rnm(ii),m_tnk(ii,k),m_tnk(ii,k),1.0D0,phi,rnm(ii),km2,m_tnk(ii,k),0.0D0,tmp_v,rnm(ii))
          call dgemm ('N','T',rnm(ii),rnm(ii),m_tnk(ii,k),1.0D0,tmp_v,rnm(ii),phi,rnm(ii),0.0D0,tmp_vcva,rnm(ii))
          deallocate(tmp_v)

          mm=sum(rnm(1:ii))-rnm(ii)
          !$OMP PARALLEL DO PRIVATE(jj, kk)
          do jj=1,rnm(ii)
            do kk=1,rnm(ii)
              !pig_vcva(k,ii*ttn-ttn+jj,ii*ttn-ttn+kk)=tmp_vcva(jj,kk)
              pig_vcva(k,mm+jj,mm+kk)=tmp_vcva(jj,kk)
            end do
          end do
          !$OMP END PARALLEL DO

          pig=0;trnx=0
          do ui=1,tn
            trnx=trnx+rnm(ui)
            trny=0
            do yi=1,ui
              trny=trny+rnm(yi)
              !$OMP PARALLEL DO PRIVATE(i,j,k2) 
              do i=1,rnm(ui)
                do j=1,rnm(yi)
                  do k2=1,mn
                    pig(trnx-rnm(ui)+i,trny-rnm(yi)+j)=pig(trnx-rnm(ui)+i,trny-rnm(yi)+j)    &
&                   +mbin(k2,yidx(ui,i),yidx(yi,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(yi)+j)
                  end do
                  if (ui.ne.yi) then
                    pig(trny-rnm(yi)+j,trnx-rnm(ui)+i)=pig(trnx-rnm(ui)+i,trny-rnm(yi)+j)
                  end if
                end do
              end do
              !$OMP END PARALLEL DO
            end do
          end do

          call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
          call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
          call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig,trn,0.0D0,py_tmp3,1)
          call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

          idum=tn
          do ll=1,mn
            do kk=1,ttn
              if (ll==k .and. kk==ii) goto 50
              idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
            end do
            do kk=2,ttn
              do mm=1,kk-1
                !idum=idum+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
              end do
            end do
          end do
50        continue
 
          aim(idum+xi,idum+xi)=tmpm(1,1)
          !print*,idum+xi,tn+nb*k-nb+xi  !'tmpm:',tmpm(1,1)/2

          !for Vg x Ve (off diagonal) !******************************************
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            pig_tmp=0.0D0
            do i=1,rnm(vi)
              pig_tmp(trnv-rnm(vi)+i,trnv-rnm(vi)+i)=1.0D0
            end do
            call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
            call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
            call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
            call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

            aim(idum+xi,vi)=tmpm(1,1)
            !print*,'vg x ve:',idum+xi,tn+nb*k-nb+xi  !tmpm(1,1)/2
          end do !vi

          !for Vg(1~[i-1]) x Vg(i) within trait
          do vi=1,xi-1
            km2=0;km2(vi,vi)=1;pig_vcva=0
            !tmp_vcva=matmul(matmul(phi,km2),transpose(phi))
            allocate(tmp_v(rnm(ii),m_tnk(ii,k)))
            call dgemm ('N','N',rnm(ii),m_tnk(ii,k),m_tnk(ii,k),1.0D0,phi,rnm(ii),km2,m_tnk(ii,k),0.0D0,tmp_v,rnm(ii))
            call dgemm ('N','T',rnm(ii),rnm(ii),m_tnk(ii,k),1.0D0,tmp_v,rnm(ii),phi,rnm(ii),0.0D0,tmp_vcva,rnm(ii))
            deallocate(tmp_v)

            mm=sum(rnm(1:ii))-rnm(ii)
            !$OMP PARALLEL DO PRIVATE (jj,kk)
            do jj=1,rnm(ii)
              do kk=1,rnm(ii)
                pig_vcva(k,mm+jj,mm+kk)=tmp_vcva(jj,kk)
              end do
            end do
            !$OMP END PARALLEL DO

            pig_tmp=0;trnx=0
            do ui=1,tn
              trnx=trnx+rnm(ui)
              trny=0
              do yi=1,ui
                trny=trny+rnm(yi)
                !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                do i=1,rnm(ui)
                  do j=1,rnm(yi)
                    do k2=1,mn
                      pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j) &
&                     +mbin(k2,yidx(ui,i),yidx(yi,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(yi)+j)
                    end do
                    if (ui.ne.yi) then
                      pig_tmp(trny-rnm(yi)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j)
                    end if
                  end do
                end do
                !$OMP END PARALLEL DO 
              end do
            end do

            call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
            call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
            call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
            call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)
            aim(idum+xi,idum+vi)=tmpm(1,1)
            !print*,'vg x vg:',idum+xi,idum+vi,tn+nb*k-nb+xi,tn+nb*k-nb+vi    !tmpm(1,1)/2
          end do


          !off diagonal between diff random eff. for Vk  !******************************
          nb3=tn
          do wi=1,k
            oo=ttn                   !for other random effect
            if (wi==k) oo=ii-1       !for within random effect
            do ii2=1,oo

              allocate(phi_2(rnm(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(ii2,wi)))
              allocate(tmp_wk_2(rnm(ii2)),tmp_vcva_2(rnm(ii2),rnm(ii2)))

              tmp_wk_2=m_wk(wi,ii2,1:rnm(ii2))
              call pol (m_tnk(ii2,wi),tmp_wk_2,rnm(ii2),phi_2)

              do vi=1,m_tnk(ii2,wi)   !check
                nb2=nb2+1

                km2_2=0;km2_2(vi,vi)=1;pig_vcva=0
                !tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
                allocate(tmp_v(rnm(ii2),m_tnk(ii2,wi)))
                call dgemm ('N','N',rnm(ii2),m_tnk(ii2,wi),m_tnk(ii2,wi),1.0D0,phi_2,rnm(ii2),km2_2,m_tnk(ii2,wi),0.0D0,tmp_v,rnm(ii2))
                call dgemm ('N','T',rnm(ii2),rnm(ii2),m_tnk(ii2,wi),1.0D0,tmp_v,rnm(ii2),phi_2,rnm(ii2),0.0D0,tmp_vcva_2,rnm(ii2))
                deallocate(tmp_v)

                mm=sum(rnm(1:ii2))-rnm(ii2)
                !$OMP PARALLEL DO PRIVATE(jj,kk) 
                do jj=1,rnm(ii2)
                  do kk=1,rnm(ii2)
                    pig_vcva(wi,mm+jj,mm+kk)=tmp_vcva_2(jj,kk)
                  end do
                end do
                !$OMP END PARALLEL DO 

                pig_tmp=0;trnx=0
                do ui=1,tn
                  trnx=trnx+rnm(ui)
                  trny=0
                  do yi=1,ui
                    trny=trny+rnm(yi)
                    !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                    do i=1,rnm(ui)
                      do j=1,rnm(yi)
                        do k2=1,mn
                          pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j)   &
&                         +mbin(k2,yidx(ui,i),yidx(yi,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(yi)+j)
                        end do
                        if (ui.ne.yi) then
                          pig_tmp(trny-rnm(yi)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j)
                        end if
                      end do
                    end do
                    !$OMP END PARALLEL DO 
                  end do
                end do

                call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
                call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
                call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
                call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

                idum2=tn
                do ll=1,mn
                  do kk=1,ttn
                    if (ll==wi .and. kk==ii2) goto 52
                    idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                  end do
                  do kk=2,ttn
                    do mm=1,kk-1
                      !idum2=idum2+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                      idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                    end do
                  end do
                end do
52              continue

                aim(idum+xi,idum2+vi)=tmpm(1,1)   !with the xi th Vg for first trait
                !print*,'*',idum+xi,idum2+vi,tn+nb*k-nb+xi,tn+nb*wi-nb+vi   !v2
              end do !vi

              nbc=0  !for no. block for covariance, i.e. 1: (2,1), 2: (3,1), 3:(3,2), ... 
              do vi2=2,m_tnk(ii2,wi)
                do ui2=1,(vi2-1)
                  nbc=nbc+1

                  km2_2=0;km2_2(ui2,vi2)=1;km2_2(vi2,ui2)=1;pig_vcva=0
                  !tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
                  allocate(tmp_v(rnm(ii2),m_tnk(ii2,wi)))
                  call dgemm ('N','N',rnm(ii2),m_tnk(ii2,wi),m_tnk(ii2,wi),1.0D0,phi_2,rnm(ii2),km2_2,m_tnk(ii2,wi),0.0D0,tmp_v,rnm(ii2))
                  call dgemm ('N','T',rnm(ii2),rnm(ii2),m_tnk(ii2,wi),1.0D0,tmp_v,rnm(ii2),phi_2,rnm(ii2),0.0D0,tmp_vcva_2,rnm(ii2))
                  deallocate(tmp_v)

                  mm=sum(rnm(1:ii2))-rnm(ii2)
                  !$OMP PARALLEL DO PRIVATE(jj,kk)
                  do jj=1,rnm(ii2)
                    do kk=1,rnm(ii2)
                      pig_vcva(wi,mm+jj,mm+kk)=tmp_vcva_2(jj,kk)
                    end do
                  end do
                  !$OMP END PARALLEL DO

                  pig_tmp=0;trnx=0
                  do ui=1,tn
                    trnx=trnx+rnm(ui)
                    trny=0
                    do vi=1,ui
                      trny=trny+rnm(vi)
                      !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                      do i=1,rnm(ui)
                        do j=1,rnm(vi)
                          do k2=1,mn
                            pig_tmp(trnx-rnm(ui)+i,trny-rnm(vi)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(vi)+j)    &
&                           +mbin(k2,yidx(ui,i),yidx(vi,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(vi)+j)
                          end do
                          if (ui.ne.vi) then
                            pig_tmp(trny-rnm(vi)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(vi)+j)
                          end if
                        end do
                      end do
                      !$OMP END PARALLEL DO 
                    end do
                  end do

                  call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
                  call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
                  call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
                  call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

                  aim(idum+xi,idum2+m_tnk(ii2,wi)+nbc)=tmpm(1,1)  !check
                  !print*,'**',idum+xi,idum2+tnk(wi)+nbc,tn+nb*k-nb+xi,tn+nb*wi-nb+tnk(wi)+nbc
                end do !ui2
              end do !vi2
              deallocate(phi_2,km2_2,tmp_wk_2,tmp_vcva_2)
            end do ! ii2

            !for Vk x CVk   *************************************
            if (wi<k) then

              do ii2=2,ttn
                do jj2=1,ii2-1

                  allocate(phi_1(rnm(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(jj2,wi)))
                  allocate(tmp_wk_1(rnm(ii2)),tmp_vcva_2(rnm(ii2),rnm(jj2)))

                  tmp_wk_1=m_wk(wi,ii2,1:rnm(ii2))
                  call pol (m_tnk(ii2,wi),tmp_wk_1,rnm(ii2),phi_1)

                  allocate(phi_2(rnm(jj2),m_tnk(jj2,wi)),tmp_wk_2(rnm(jj2)))
                  tmp_wk_2=m_wk(wi,jj2,1:rnm(jj2))
                  call pol (m_tnk(jj2,wi),tmp_wk_2,rnm(jj2),phi_2)

                  nbc2=0
                  do vi=1,m_tnk(ii2,wi)
                    do vi2=1,m_tnk(jj2,wi)
                      nbc2=nbc2+1
                      km2_2=0;km2_2(vi,vi2)=1;pig_vcva=0
                      !tmp_vcva_2=matmul(matmul(phi_1,km2_2),transpose(phi_2))
                      allocate(tmp_v(rnm(ii2),m_tnk(jj2,wi)))
                      call dgemm ('N','N',rnm(ii2),m_tnk(jj2,wi),m_tnk(ii2,wi),1.0D0,phi_1,rnm(ii2),km2_2,m_tnk(ii2,wi),0.0D0,tmp_v,rnm(ii2))
                      call dgemm ('N','T',rnm(ii2),rnm(jj2),m_tnk(jj2,wi),1.0D0,tmp_v,rnm(ii2),phi_2,rnm(jj2),0.0D0,tmp_vcva_2,rnm(ii2))
                      deallocate(tmp_v)

                      mm=sum(rnm(1:ii2))-rnm(ii2)
                      mm2=sum(rnm(1:jj2))-rnm(jj2)
                      !$OMP PARALLEL DO PRIVATE(ll,kk)
                      do ll=1,rnm(ii2)
                        do kk=1,rnm(jj2)
                          pig_vcva(wi,mm+ll,mm2+kk)=tmp_vcva_2(ll,kk)
                          pig_vcva(wi,mm2+kk,mm+ll)=tmp_vcva_2(ll,kk)
                        end do
                      end do
                      !$OMP END PARALLEL DO

                      pig_tmp=0;trnx=0
                      do ui=1,tn
                        trnx=trnx+rnm(ui)
                        trny=0
                        do yi=1,ui
                          trny=trny+rnm(yi)
                          !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                          do i=1,rnm(ui)
                            do j=1,rnm(yi)
                              do k2=1,mn
                                pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j)&
&                               +mbin(k2,yidx(ui,i),yidx(yi,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(yi)+j)
                              end do
                              if (ui.ne.yi) then
                                pig_tmp(trny-rnm(yi)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j)
                              end if
                            end do
                          end do
                          !$OMP END PARALLEL DO 
                        end do
                      end do
                      call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
                      call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
                      call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
                      call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

                      idum2=tn
                      do ll=1,mn
                        do kk=1,ttn
                          idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                        end do
                        do kk=2,ttn
                          do mm=1,kk-1
                            if (ll==wi .and. kk==ii2 .and. mm==jj2) goto 60
                            !idum2=idum2+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                            idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                          end do
                        end do
                      end do
60                    continue

                      aim(idum+xi,idum2+nbc2)=tmpm(1,1)  !check
                      !print*,'Vk x CVk in prv',idum+xi,idum2+nbc2
                    end do ! vi2
                  end do ! vi

                  deallocate(phi_1,phi_2,km2_2)
                  deallocate(tmp_wk_1,tmp_wk_2,tmp_vcva_2)

                end do ! jj2
              end do ! ii2
            end if ! wk<k

          end do !wi

          !dldv Vg *********************************************************

          tr1=0
          do i=1,trn
            do j=1,trn
              tr1=tr1+pm(i,j)*pig(j,i)
            end do
          end do
          tr1=-0.5*tr1

          !tmpm=matmul(matmul(transpose(py),pig),py)
          call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
          call dgemm ('N','N',1,1,trn,1.0D0,py_tmp1,1,py,trn,0.0D0,tmpm,1)
        
          dldv(idum+xi,1)=tmpm(1,1)*0.5+tr1 !************************************
          !print*,'dldv:',tnk+nb*k-nb+xi,1,dldv(tn+nb*k-nb+xi,1),tn+nb*k-nb+xi
        end do !xi

        !for cov (diagonal) *************************************************
        nbc=0  !for no. block for covariance, i.e. 1: (2,1), 2: (3,1), 3:(3,2),... 
        do xi=2,m_tnk(ii,k)
          do yi=1,xi-1
            nbc=nbc+1
            km2=0;km2(xi,yi)=1;km2(yi,xi)=1;pig_vcva=0
            !tmp_vcva=matmul(matmul(phi,km2),transpose(phi))

            allocate(tmp_v(rnm(ii),m_tnk(ii,k)))
            call dgemm ('N','N',rnm(ii),m_tnk(ii,k),m_tnk(ii,k),1.0D0,phi,rnm(ii),km2,m_tnk(ii,k),0.0D0,tmp_v,rnm(ii))
            call dgemm ('N','T',rnm(ii),rnm(ii),m_tnk(ii,k),1.0D0,tmp_v,rnm(ii),phi,rnm(ii),0.0D0,tmp_vcva,rnm(ii))
            deallocate(tmp_v)

            mm=sum(rnm(1:ii))-rnm(ii)
            !$OMP PARALLEL DO PRIVATE(jj,kk)
            do jj=1,rnm(ii)
              do kk=1,rnm(ii)
                pig_vcva(k,mm+jj,mm+kk)=tmp_vcva(jj,kk)
              end do
            end do
            !$OMP END PARALLEL DO

          pig=0;trnx=0
          do ui=1,tn
            trnx=trnx+rnm(ui)
            trny=0
            do vi=1,ui
              trny=trny+rnm(vi)
              !$OMP PARALLEL DO PRIVATE(i,j,k2) 
              do i=1,rnm(ui)
                do j=1,rnm(vi)
                  do k2=1,mn
                    pig(trnx-rnm(ui)+i,trny-rnm(vi)+j)=pig(trnx-rnm(ui)+i,trny-rnm(vi)+j)     &
&                   +mbin(k2,yidx(ui,i),yidx(vi,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(vi)+j)
                  end do
                  if (ui.ne.vi) then
                    pig(trny-rnm(vi)+j,trnx-rnm(ui)+i)=pig(trnx-rnm(ui)+i,trny-rnm(vi)+j)
                  end if
                end do
              end do
              !$OMP END PARALLEL DO 
            end do
          end do

          call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
          call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
          call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig,trn,0.0D0,py_tmp3,1)
          call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

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
            aim(idum+m_tnk(ii,k)+nbc,idum+m_tnk(ii,k)+nbc)=tmpm(1,1)     !check
            !print*,'cov***:',idum+tnk(k)+nbc,tn+nb*k-nb+tnk(k)+nbc

            !for cov x Ve (off diagonal) ************************************
            trnv=0
            do vi=1,tn
              trnv=trnv+rnm(vi)

              pig_tmp=0.0D0
              do i=1,rnm(vi)
                pig_tmp(trnv-rnm(vi)+i,trnv-rnm(vi)+i)=1.0D0
              end do
              call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
              call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
              call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
              call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

              aim(idum+m_tnk(ii,k)+nbc,vi)=tmpm(1,1)          !check with Ve1
              !print*,'cov x
            end do

            !for cov x Vg **************************************************
            do vi=1,m_tnk(ii,k)

              km2=0;km2(vi,vi)=1;pig_vcva=0
              !tmp_vcva=matmul(matmul(phi,km2),transpose(phi))

              allocate(tmp_v(rnm(ii),m_tnk(ii,k)))
              call dgemm ('N','N',rnm(ii),m_tnk(ii,k),m_tnk(ii,k),1.0D0,phi,rnm(ii),km2,m_tnk(ii,k),0.0D0,tmp_v,rnm(ii))
              call dgemm ('N','T',rnm(ii),rnm(ii),m_tnk(ii,k),1.0D0,tmp_v,rnm(ii),phi,rnm(ii),0.0D0,tmp_vcva,rnm(ii))
              deallocate(tmp_v)

              mm=sum(rnm(1:ii))-rnm(ii)
             !$OMP PARALLEL DO PRIVATE(jj,kk)
              do jj=1,rnm(ii)
                do kk=1,rnm(ii)
                  pig_vcva(k,mm+jj,mm+kk)=tmp_vcva(jj,kk)
                end do
              end do
             !$OMP END PARALLEL DO

              pig_tmp=0;trnx=0
              do ui=1,tn
                trnx=trnx+rnm(ui)
                trny=0
                do wi=1,ui
                  trny=trny+rnm(wi)
                  !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                  do i=1,rnm(ui)
                    do j=1,rnm(wi)
                      do k2=1,mn
                        pig_tmp(trnx-rnm(ui)+i,trny-rnm(wi)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(wi)+j)  &
&                       +mbin(k2,yidx(ui,i),yidx(wi,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(wi)+j)
                      end do
                      if (ui.ne.wi) then
                        pig_tmp(trny-rnm(wi)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(wi)+j)
                      end if
                    end do
                  end do
              !$OMP END PARALLEL DO
                end do
              end do
              call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
              call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
              call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
              call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

              aim(idum+m_tnk(ii,k)+nbc,idum+vi)=tmpm(1,1)  !with the same Vg1 for first trait
              !print*,'cov x vg:',idum+tnk(k)+nbc,idum+vi,tn+nb*k-nb+tnk(k)+nbc,tn+nb*k-nb+vi
            end do

            !for cov x cov *************************************************
            trnv=rnm(1)
            nbc2=0    !no block for cov (vi,ui)
            do vi=2,m_tnk(ii,k)
              trnv=trnv+rnm(vi)
              trnu=0
              do ui=1,(vi-1)
                nbc2=nbc2+1

                if (nbc2<nbc) then  !off diagonal within cova ************************
                  trnu=trnu+rnm(ui)

                  km2=0;km2(vi,ui)=1;km2(ui,vi)=1;pig_vcva=0
                  !tmp_vcva=matmul(matmul(phi,km2),transpose(phi))
                  allocate(tmp_v(rnm(ii),m_tnk(ii,k)))
                  call dgemm ('N','N',rnm(ii),m_tnk(ii,k),m_tnk(ii,k),1.0D0,phi,rnm(ii),km2,m_tnk(ii,k),0.0D0,tmp_v,rnm(ii))
                  call dgemm ('N','T',rnm(ii),rnm(ii),m_tnk(ii,k),1.0D0,tmp_v,rnm(ii),phi,rnm(ii),0.0D0,tmp_vcva,rnm(ii))
                  deallocate(tmp_v)

                  mm=sum(rnm(1:ii))-rnm(ii)
                  !$OMP PARALLEL DO PRIVATE(jj,kk)
                  do jj=1,rnm(ii)
                    do kk=1,rnm(ii)
                      pig_vcva(k,mm+jj,mm+kk)=tmp_vcva(jj,kk)
                    end do
                  end do
                  !$OMP END PARALLEL DO

                  pig_tmp=0;trnx=0
                  do wj=1,tn
                    trnx=trnx+rnm(wj)
                    trny=0
                    do wi=1,wj
                      trny=trny+rnm(wi)
                      !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                      do i=1,rnm(wj)
                        do j=1,rnm(wi)
                          do k2=1,mn
                            pig_tmp(trnx-rnm(wj)+i,trny-rnm(wi)+j)=pig_tmp(trnx-rnm(wj)+i,trny-rnm(wi)+j) &
&                           +mbin(k2,yidx(wj,i),yidx(wi,j))*pig_vcva(k2,trnx-rnm(wj)+i,trny-rnm(wi)+j)
                          end do
                          if (wj.ne.wi) then
                            pig_tmp(trny-rnm(wi)+j,trnx-rnm(wj)+i)=pig_tmp(trnx-rnm(wj)+i,trny-rnm(wi)+j)
                          end if
                        end do
                      end do
                      !$OMP END PARALLEL DO 
                    end do
                  end do
                  call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
                  call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
                  call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
                  call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

                  aim(idum+m_tnk(ii,k)+nbc,idum+m_tnk(ii,k)+nbc2)=tmpm(1,1) !with the wi th cov for first trait
                  !print*,'cov x cov:',idum+tnk(k)+nbc,idum+tnk(k)+nbc2,tn+nb*k-nb+tnk(k)+nbc,tn+nb*k-nb+tnk(k)+nbc2
                end if ! *************************************************************
              end do !ui
            end do !vi

            !off diagonal berween other random eff. for cov   !****************************
            do wi=1,k
              oo=ttn                   !for other random effect
              if (wi==k) oo=ii-1       !for within random effect

              do ii2=1,oo
                allocate(phi_2(rnm(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(ii2,wi)))
                allocate(tmp_wk_2(rnm(ii2)),tmp_vcva_2(rnm(ii2),rnm(ii2)))

                tmp_wk_2=m_wk(wi,ii2,1:rnm(ii2))
                call pol (m_tnk(ii2,wi),tmp_wk_2,rnm(ii2),phi_2)

                do vi=1,m_tnk(ii2,wi)   !check
                  nb2=nb2+1

                  km2_2=0;km2_2(vi,vi)=1;pig_vcva=0
                  !tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))

                  allocate(tmp_v(rnm(ii2),m_tnk(ii2,wi)))
                  call dgemm ('N','N',rnm(ii2),m_tnk(ii2,wi),m_tnk(ii2,wi),1.0D0,phi_2,rnm(ii2),km2_2,m_tnk(ii2,wi),0.0D0,tmp_v,rnm(ii2))
                  call dgemm ('N','T',rnm(ii2),rnm(ii2),m_tnk(ii2,wi),1.0D0,tmp_v,rnm(ii2),phi_2,rnm(ii2),0.0D0,tmp_vcva_2,rnm(ii2))
                  deallocate(tmp_v)

                  mm=sum(rnm(1:ii2))-rnm(ii2)
                  !$OMP PARALLEL DO PRIVATE(jj,kk)
                  do jj=1,rnm(ii2)
                    do kk=1,rnm(ii2)
                      pig_vcva(wi,mm+jj,mm+kk)=tmp_vcva_2(jj,kk)
                    end do
                  end do
                  !$OMP END PARALLEL DO 

              pig_tmp=0;trnx=0
              do ui=1,tn
                trnx=trnx+rnm(ui)
                trny=0
                do yi2=1,ui
                  trny=trny+rnm(yi2)
                  !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                  do i=1,rnm(ui)
                    do j=1,rnm(yi2)
                      do k2=1,mn
                        pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j)      &
&                       +mbin(k2,yidx(ui,i),yidx(yi2,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(yi2)+j)
                      end do
                      if (ui.ne.yi2) then
                        pig_tmp(trny-rnm(yi2)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j)
                      end if
                    end do
                  end do
                  !$OMP END PARALLEL DO
                end do
              end do

              call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
              call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
              call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
              call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

                idum2=tn
                do ll=1,mn
                  do kk=1,ttn
                    if (ll==wi .and. kk==ii2) goto 54
                    idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                  end do
                  do kk=2,ttn
                    do mm=1,kk-1
                      idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                    end do
                  end do
                end do
54              continue

                aim(idum+m_tnk(ii,k)+nbc,idum2+vi)=tmpm(1,1)  !with the wi th Vg for first trait
                !print*,'****',idum+tnk(k)+nbc,idum2+vi,tn+nb*k-nb+tnk+nbc,tn+nb*wi-nb+vi
              end do

              nbc2=0    !no block for cov (vi,ui)
                do vi2=2,m_tnk(ii2,wi)
                  do ui2=1,(vi2-1)
                    nbc2=nbc2+1

                    km2_2=0;km2_2(ui2,vi2)=1;km2_2(vi2,ui2)=1;pig_vcva=0
                    !tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))

                    allocate(tmp_v(rnm(ii2),m_tnk(ii2,wi)))
                    call dgemm ('N','N',rnm(ii2),m_tnk(ii2,wi),m_tnk(ii2,wi),1.0D0,phi_2,rnm(ii2),km2_2,m_tnk(ii2,wi),0.0D0,tmp_v,rnm(ii2))
                    call dgemm ('N','T',rnm(ii2),rnm(ii2),m_tnk(ii2,wi),1.0D0,tmp_v,rnm(ii2),phi_2,rnm(ii2),0.0D0,tmp_vcva_2,rnm(ii2))
                    deallocate(tmp_v)

                    mm=sum(rnm(1:ii2))-rnm(ii2)
                    !$OMP PARALLEL DO PRIVATE(jj,kk)
                    do jj=1,rnm(ii2)
                      do kk=1,rnm(ii2)
                        pig_vcva(wi,mm+jj,mm+kk)=tmp_vcva_2(jj,kk)
                      end do
                    end do
                    !$OMP END PARALLEL DO 

                pig_tmp=0;trnx=0
                do ui=1,tn
                  trnx=trnx+rnm(ui)
                  trny=0
                  do vi=1,ui
                    trny=trny+rnm(vi)
                    !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                    do i=1,rnm(ui)
                      do j=1,rnm(vi)
                        do k2=1,mn
                          pig_tmp(trnx-rnm(ui)+i,trny-rnm(vi)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(vi)+j)    &
&                         +mbin(k2,yidx(ui,i),yidx(vi,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(vi)+j)
                        end do
                        if (ui.ne.vi) then
                          pig_tmp(trny-rnm(vi)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(vi)+j)
                        end if
                      end do
                    end do
                    !$OMP END PARALLEL DO
                  end do
                end do

                call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
                call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
                call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
                call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

                  aim(idum+m_tnk(ii,k)+nbc,idum2+m_tnk(ii2,wi)+nbc2)=tmpm(1,1) !with the wi th cov for first trait

                  !print*,'*****',tn+nb*k-nb+tnk+nbc,tn+nb*wi-nb+tnk+nbc2
                end do !ui2
              end do !vi2

              deallocate(phi_2,km2_2,tmp_wk_2,tmp_vcva_2)
            end do ! ii2

            !for cov x CVk   *************************************
            if (wi<k) then

              do ii2=2,ttn
                do jj2=1,ii2-1

                  allocate(phi_1(rnm(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(jj2,wi)))
                  allocate(tmp_wk_1(rnm(ii2)),tmp_vcva_2(rnm(ii2),rnm(jj2)))

                  tmp_wk_1=m_wk(wi,ii2,1:rnm(ii2))
                  call pol (m_tnk(ii2,wi),tmp_wk_1,rnm(ii2),phi_1)

                  allocate(phi_2(rnm(jj2),m_tnk(jj2,wi)),tmp_wk_2(rnm(jj2)))
                  tmp_wk_2=m_wk(wi,jj2,1:rnm(jj2))
                  call pol (m_tnk(jj2,wi),tmp_wk_2,rnm(jj2),phi_2)

                  nbc2=0
                  do vi=1,m_tnk(ii2,wi)
                    do vi2=1,m_tnk(jj2,wi)
                      nbc2=nbc2+1

                    km2_2=0;km2_2(vi,vi2)=1;pig_vcva=0
                    !tmp_vcva_2=matmul(matmul(phi_1,km2_2),transpose(phi_2))

                    allocate(tmp_v(rnm(ii2),m_tnk(jj2,wi)))
                    call dgemm ('N','N',rnm(ii2),m_tnk(jj2,wi),m_tnk(ii2,wi),1.0D0,phi_1,rnm(ii2),km2_2,m_tnk(ii2,wi),0.0D0,tmp_v,rnm(ii2))
                    call dgemm ('N','T',rnm(ii2),rnm(jj2),m_tnk(jj2,wi),1.0D0,tmp_v,rnm(ii2),phi_2,rnm(jj2),0.0D0,tmp_vcva_2,rnm(ii2))
                    deallocate(tmp_v)

                    mm=sum(rnm(1:ii2))-rnm(ii2)
                    mm2=sum(rnm(1:jj2))-rnm(jj2)
                    !$OMP PARALLEL DO PRIVATE(ll,kk) 
                    do ll=1,rnm(ii2)
                      do kk=1,rnm(jj2)
                        pig_vcva(wi,mm+ll,mm2+kk)=tmp_vcva_2(ll,kk)
                        pig_vcva(wi,mm2+kk,mm+ll)=tmp_vcva_2(ll,kk)
                      end do
                    end do
                    !$OMP END PARALLEL DO

                pig_tmp=0;trnx=0
                do ui=1,tn
                  trnx=trnx+rnm(ui)
                  trny=0
                  do yi2=1,ui
                    trny=trny+rnm(yi2)
                    !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                    do i=1,rnm(ui)
                      do j=1,rnm(yi2)
                        do k2=1,mn
                          pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j) &
&                         +mbin(k2,yidx(ui,i),yidx(yi2,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(yi2)+j)
                        end do
                        if (ui.ne.yi2) then
                          pig_tmp(trny-rnm(yi2)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j)
                        end if
                      end do
                    end do
                  !$OMP END PARALLEL DO
                  end do
                end do

                call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
                call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
                call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
                call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

                    idum2=tn
                    do ll=1,mn
                      do kk=1,ttn
                        idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                      end do
                      do kk=2,ttn
                        do mm=1,kk-1
                          if (ll==wi .and. kk==ii2 .and. mm==jj2) goto 61
                          !idum2=idum2+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                          idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                        end do
                      end do
                    end do
61                  continue

                    aim(idum+m_tnk(ii,k)+nbc,idum2+nbc2)=tmpm(1,1)
                    !print*,'cov x CVk in prv',idum+m_tnk(ii,k)+nbc,idum2+nbc2

                    end do !vi2
                  end do ! vi

                  deallocate(phi_1,phi_2,km2_2)
                  deallocate(tmp_wk_1,tmp_wk_2,tmp_vcva_2)
                  !print*,'??'

                end do ! jj2
              end do ! ii2
            end if ! wk<k
                                                                                      
          end do !wi

          ! dldv for cov
          tr1=0
          do i=1,trn
            do j=1,trn
              tr1=tr1+pm(i,j)*pig(j,i)
            end do
          end do
          tr1=-0.5*tr1

          call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
          call dgemm ('N','N',1,1,trn,1.0D0,py_tmp1,1,py,trn,0.0D0,tmpm,1)

          dldv(idum+m_tnk(ii,k)+nbc,1)=tmpm(1,1)*0.5+tr1
          !print*,'dldv:',dldv(tn+nb*k-nb+tnk+nbc,1),tn+nb*k-nb+tnk+nbc
        end do !yi
      end do !xi

      deallocate(phi,km2,tmp_wk,tmp_vcva)

    end do !ii


      !CVk, i.e. covariance between mrrm traits  **************************
      do ii=2,ttn
        do jj=1,ii-1
          allocate(phi(rnm(ii),m_tnk(ii,k)),km2(m_tnk(ii,k),m_tnk(jj,k)))
          allocate(tmp_wk(rnm(ii)),tmp_vcva(rnm(ii),rnm(jj)))

          tmp_wk=m_wk(k,ii,1:rnm(ii))
          call pol (m_tnk(ii,k),tmp_wk,rnm(ii),phi)

          allocate(phi2(rnm(jj),m_tnk(jj,k)),tmp_wk2(rnm(jj)))
          tmp_wk2=m_wk(k,jj,1:rnm(jj))
          call pol (m_tnk(jj,k),tmp_wk2,rnm(jj),phi2)

          nbc=0
          do xi=1,m_tnk(ii,k)
          do xi2=1,m_tnk(jj,k)
            nbc=nbc+1

            km2=0;km2(xi,xi2)=1;pig_vcva=0
            !tmp_vcva=matmul(matmul(phi,km2),transpose(phi2))
            allocate(tmp_v(rnm(ii),m_tnk(jj,k)))
            call dgemm ('N','N',rnm(ii),m_tnk(jj,k),m_tnk(ii,k),1.0D0,phi,rnm(ii),km2,m_tnk(ii,k),0.0D0,tmp_v,rnm(ii))
            call dgemm ('N','T',rnm(ii),rnm(jj),m_tnk(jj,k),1.0D0,tmp_v,rnm(ii),phi2,rnm(jj),0.0D0,tmp_vcva,rnm(ii))
            deallocate(tmp_v)

            mm=sum(rnm(1:ii))-rnm(ii)
            mm2=sum(rnm(1:jj))-rnm(jj)
            !$OMP PARALLEL DO PRIVATE(ll,kk)
            do ll=1,rnm(ii)
              do kk=1,rnm(jj)
                pig_vcva(k,mm+ll,mm2+kk)=tmp_vcva(ll,kk)
                pig_vcva(k,mm2+kk,mm+ll)=tmp_vcva(ll,kk)
              end do
            end do
            !$OMP END PARALLEL DO

            pig=0;trnx=0
            do ui=1,tn
              trnx=trnx+rnm(ui)
              trny=0
              do yi=1,ui
                trny=trny+rnm(yi)
                !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                do i=1,rnm(ui)
                  do j=1,rnm(yi)
                    do k2=1,mn
                      pig(trnx-rnm(ui)+i,trny-rnm(yi)+j)=pig(trnx-rnm(ui)+i,trny-rnm(yi)+j)&
&                     +mbin(k2,yidx(ui,i),yidx(yi,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(yi)+j)
                    end do
                    if (ui.ne.yi) then
                      pig(trny-rnm(yi)+j,trnx-rnm(ui)+i)=pig(trnx-rnm(ui)+i,trny-rnm(yi)+j)
                    end if
                  end do
                end do
                !$OMP END PARALLEL DO
              end do
            end do

            call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
            call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
            call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig,trn,0.0D0,py_tmp3,1)
            call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

            idum=tn
            do ll=1,mn
              do kk=1,ttn
                idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
              end do
              do kk=2,ttn
                  do mm=1,kk-1
                  if (ll==k .and. kk==ii .and. mm==jj) goto 65
                    !idum=idum+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                    idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
                  end do
                end do

              end do
65            continue

            aim(idum+nbc,idum+nbc)=tmpm(1,1)
            !print*,'CVg',idum+xi,tmpm(1,1)

            !for CVg x Ve (off diagonal)
            !!******************************************
            trnv=0
            do vi=1,tn
              trnv=trnv+rnm(vi)
              pig_tmp=0.0D0
              do i=1,rnm(vi)
                pig_tmp(trnv-rnm(vi)+i,trnv-rnm(vi)+i)=1.0D0
              end do
              call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
              call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
              call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
              call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

              aim(idum+nbc,vi)=tmpm(1,1)
              !print*,'Cvg x ve:',idum+xi,vi,tmpm(1,1)
            end do ! vi

              !for CVk(1~[i-1]) x CVk(i) *********************************
              nbc2=0
              do vi=1,m_tnk(ii,k)
              do vi2=1,m_tnk(jj,k)
                nbc2=nbc2+1
                if (nbc2<nbc) then
                km2=0;km2(vi,vi2)=1;pig_vcva=0
                !tmp_vcva=matmul(matmul(phi,km2),transpose(phi2))
                allocate(tmp_v(rnm(ii),m_tnk(jj,k)))
                call dgemm ('N','N',rnm(ii),m_tnk(jj,k),m_tnk(ii,k),1.0D0,phi,rnm(ii),km2,m_tnk(ii,k),0.0D0,tmp_v,rnm(ii))
                call dgemm ('N','T',rnm(ii),rnm(jj),m_tnk(jj,k),1.0D0,tmp_v,rnm(ii),phi2,rnm(jj),0.0D0,tmp_vcva,rnm(ii))
                deallocate(tmp_v)

                mm=sum(rnm(1:ii))-rnm(ii)
                mm2=sum(rnm(1:jj))-rnm(jj)
                !$OMP PARALLEL DO PRIVATE(ll,kk)
                do ll=1,rnm(ii)
                  do kk=1,rnm(jj)
                    pig_vcva(k,mm+ll,mm2+kk)=tmp_vcva(ll,kk)
                    pig_vcva(k,mm2+kk,mm+ll)=tmp_vcva(ll,kk)
                  end do
                end do
                !$OMP END PARALLEL DO

              pig_tmp=0;trnx=0
              do ui=1,tn
                trnx=trnx+rnm(ui)
                trny=0
                do yi=1,ui
                  trny=trny+rnm(yi)
                  !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                  do i=1,rnm(ui)
                    do j=1,rnm(yi)
                      do k2=1,mn
                        pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j) &
&                       +mbin(k2,yidx(ui,i),yidx(yi,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(yi)+j)

                      end do
                      if (ui.ne.yi) then
                        pig_tmp(trny-rnm(yi)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi)+j)
                      end if
                    end do
                  end do
                  !$OMP END PARALLEL DO
                end do
              end do

              call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
              call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
              call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
              call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)
              aim(idum+nbc,idum+nbc2)=tmpm(1,1)
              !print*,'Cvg x Cvg:',idum+vi,tmpm(1,1)
                end if
              end do !vi2
              end do !vi

            !off diagonal between diff random eff. for CVk
            !******************************
            nb3=tn
              do wi=1,k     !check in case multiple random effects
                do ii2=2,ttn
                  do jj2=1,ii2-1
                    if (wi==k .and. ii2>ii) then
                    elseif (wi==k .and. ii2==ii .and. jj2>=jj) then
                    else


                      allocate(phi_1(rnm(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(jj2,wi)))
                      allocate(tmp_wk_1(rnm(ii2)),tmp_vcva_2(rnm(ii2),rnm(jj2)))

                      tmp_wk_1=m_wk(wi,ii2,1:rnm(ii2))
                      call pol (m_tnk(ii2,wi),tmp_wk_1,rnm(ii2),phi_1)

                      allocate(phi_2(rnm(jj2),m_tnk(jj2,wi)),tmp_wk_2(rnm(jj2)))
                      tmp_wk_2=m_wk(wi,jj2,1:rnm(jj2))
                      call pol (m_tnk(jj2,wi),tmp_wk_2,rnm(jj2),phi_2)

                      nbc2=0
                      do vi=1,m_tnk(ii2,wi)
                      do vi2=1,m_tnk(jj2,wi)
                        nbc2=nbc2+1
                        km2_2=0;km2_2(vi,vi2)=1;pig_vcva=0
                        !tmp_vcva_2=matmul(matmul(phi_1,km2_2),transpose(phi_2))
                        allocate(tmp_v(rnm(ii2),m_tnk(jj2,wi)))
                        call dgemm ('N','N',rnm(ii2),m_tnk(jj2,wi),m_tnk(ii2,wi),1.0D0,phi_1,rnm(ii2),km2_2,m_tnk(ii2,wi),0.0D0,tmp_v,rnm(ii2))
                        call dgemm ('N','T',rnm(ii2),rnm(jj2),m_tnk(jj2,wi),1.0D0,tmp_v,rnm(ii2),phi_2,rnm(jj2),0.0D0,tmp_vcva_2,rnm(ii2))
                        deallocate(tmp_v)

                        mm=sum(rnm(1:ii2))-rnm(ii2)
                        mm2=sum(rnm(1:jj2))-rnm(jj2)
                        !$OMP PARALLEL DO PRIVATE(ll,kk)
                        do ll=1,rnm(ii2)
                          do kk=1,rnm(jj2)
                            pig_vcva(wi,mm+ll,mm2+kk)=tmp_vcva_2(ll,kk)
                            pig_vcva(wi,mm2+kk,mm+ll)=tmp_vcva_2(ll,kk)
                          end do
                        end do
                        !$OMP END PARALLEL DO

                pig_tmp=0;trnx=0
                do ui=1,tn
                  trnx=trnx+rnm(ui)
                  trny=0
                  do yi2=1,ui
                    trny=trny+rnm(yi2)
                    !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                    do i=1,rnm(ui)
                      do j=1,rnm(yi2)
                        do k2=1,mn
                          pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j)&
&                         +mbin(k2,yidx(ui,i),yidx(yi2,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(yi2)+j)
                        end do
                        if (ui.ne.yi2) then
                          pig_tmp(trny-rnm(yi2)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j)
                        end if
                      end do
                    end do
                  !$OMP END PARALLEL DO
                  end do
                end do

                call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
                call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
                call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
                call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

                        idum2=tn
                        do ll=1,mn
                          do kk=1,ttn
                            idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                          end do
                          do kk=2,ttn
                            do mm=1,kk-1
                              if (ll==wi .and. kk==ii2 .and. mm==jj2) goto 70
                              !idum2=idum2+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                              idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                            end do
                          end do
                        end do
70                      continue

                        !aim(idum+xi,idum2+vi)=v2
                        aim(idum+nbc,idum2+nbc2)=tmpm(1,1)
                        !print*,'CVk x CVk in
                        !prv',idum+nbc,idum2+nbc2,idum2,nbc2,v2
                      end do  ! vi2
                      end do ! vi

                      deallocate (phi_1,phi_2,km2_2,tmp_wk_1,tmp_vcva_2,tmp_wk_2)

                    end if
                  end do ! jj2
                end do ! ii2
              end do  ! wi

              ! CVk x Vk, off-diagonal between other traits 
              do wi=1,k
                do ii2=1,ttn
                  allocate(phi_2(rnm(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(ii2,wi)))
                  allocate(tmp_wk_2(rnm(ii2)),tmp_vcva_2(rnm(ii2),rnm(ii2)))

                  tmp_wk_2=m_wk(wi,ii2,1:rnm(ii2))
                  call pol (m_tnk(ii2,wi),tmp_wk_2,rnm(ii2),phi_2)

                  do vi=1,m_tnk(ii2,wi)   !check
                   nb2=nb2+1
                   km2_2=0;km2_2(vi,vi)=1;pig_vcva=0
                   !tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))

                    allocate(tmp_v(rnm(ii2),m_tnk(ii2,wi)))
                    call dgemm ('N','N',rnm(ii2),m_tnk(ii2,wi),m_tnk(ii2,wi),1.0D0,phi_2,rnm(ii2),km2_2,m_tnk(ii2,wi),0.0D0,tmp_v,rnm(ii2))
                    call dgemm ('N','T',rnm(ii2),rnm(ii2),m_tnk(ii2,wi),1.0D0,tmp_v,rnm(ii2),phi_2,rnm(ii2),0.0D0,tmp_vcva_2,rnm(ii2))
                    deallocate(tmp_v)

                    mm=sum(rnm(1:ii2))-rnm(ii2)
                    !$OMP PARALLEL DO PRIVATE(jj2,kk) 
                    do jj2=1,rnm(ii2)
                      do kk=1,rnm(ii2)
                        pig_vcva(wi,mm+jj2,mm+kk)=tmp_vcva_2(jj2,kk)
                      end do
                    end do
                    !$OMP END PARALLEL DO

                  pig_tmp=0;trnx=0
                  do ui=1,tn
                    trnx=trnx+rnm(ui)
                    trny=0

                    do vi2=1,ui
                      trny=trny+rnm(vi2)
                      !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                      do i=1,rnm(ui)
                        do j=1,rnm(vi2)
                          do k2=1,mn
                            pig_tmp(trnx-rnm(ui)+i,trny-rnm(vi2)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(vi2)+j) &
&                           +mbin(k2,yidx(ui,i),yidx(vi2,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(vi2)+j)
                          end do
                          if (ui.ne.vi2) then
                            pig_tmp(trny-rnm(vi2)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(vi2)+j)
                          end if
                        end do
                      end do
                      !$OMP END PARALLEL DO
                    end do
                  end do

                  call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
                  call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
                  call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
                  call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

                    idum2=tn
                    do ll=1,mn
                      do kk=1,ttn
                        if (ll==wi .and. kk==ii2) goto 58
                        idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                      end do
                      do kk=2,ttn
                        do mm=1,kk-1
                          !idum2=idum2+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                          idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                        end do
                      end do
                    end do
58                  continue

                    aim(idum+nbc,idum2+vi)=tmpm(1,1)
                    !print*,'CVk x Vk',idum+nbc,idum2+vi,idum2,vi,v2
                  end do ! vi

                  nbc2=0  ! CVk x cov, i.e. 1: (2,1),2:(3,1),3:(3,2),... 
                  do vi2=2,m_tnk(ii2,wi)
                    do ui2=1,(vi2-1)
                      nbc2=nbc2+1
                      km2_2=0;km2_2(ui2,vi2)=1;km2_2(vi2,ui2)=1;pig_vcva=0
                      !tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))

                      allocate(tmp_v(rnm(ii2),m_tnk(ii2,wi)))
                      call dgemm ('N','N',rnm(ii2),m_tnk(ii2,wi),m_tnk(ii2,wi),1.0D0,phi_2,rnm(ii2),km2_2,m_tnk(ii2,wi),0.0D0,tmp_v,rnm(ii2))
                      call dgemm ('N','T',rnm(ii2),rnm(ii2),m_tnk(ii2,wi),1.0D0,tmp_v,rnm(ii2),phi_2,rnm(ii2),0.0D0,tmp_vcva_2,rnm(ii2))
                      deallocate(tmp_v)

                      mm=sum(rnm(1:ii2))-rnm(ii2)
                      !$OMP PARALLEL DO PRIVATE(jj2,kk)
                      do jj2=1,rnm(ii2)
                        do kk=1,rnm(ii2)
                          pig_vcva(wi,mm+jj2,mm+kk)=tmp_vcva_2(jj2,kk)
                        end do
                      end do
                      !$OMP END PARALLEL DO

                pig_tmp=0;trnx=0
                do ui=1,tn
                  trnx=trnx+rnm(ui)
                  trny=0
                  do yi2=1,ui
                    trny=trny+rnm(yi2)
                    !$OMP PARALLEL DO PRIVATE(i,j,k2) 
                    do i=1,rnm(ui)
                      do j=1,rnm(yi2)
                        do k2=1,mn
                          pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j) &
&                         +mbin(k2,yidx(ui,i),yidx(yi2,j))*pig_vcva(k2,trnx-rnm(ui)+i,trny-rnm(yi2)+j)
                        end do
                        if (ui.ne.yi2) then
                          pig_tmp(trny-rnm(yi2)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j)

                          pig_tmp(trny-rnm(yi2)+j,trnx-rnm(ui)+i)=pig_tmp(trnx-rnm(ui)+i,trny-rnm(yi2)+j)
                        end if
                      end do
                    end do
                  !$OMP END PARALLEL DO
                  end do
                end do

                call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
                call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp1,1,pm,trn,0.0D0,py_tmp2,1)
                call dgemm ('N','N',1,trn,trn,1.0D0,py_tmp2,1,pig_tmp,trn,0.0D0,py_tmp3,1)
                call dgemm ('N','N',1,1,trn,1.0D0,py_tmp3,1,py,trn,0.0D0,tmpm,1)

                    aim(idum+nbc,idum2+m_tnk(ii2,wi)+nbc2)=tmpm(1,1)
                    !print*,'CVk x cov',idum+nbc,idum2+m_tnk(ii2,wi)+nbc2
                  end do !ui2
                end do !vi2

                deallocate(phi_2,km2_2,tmp_wk_2,tmp_vcva_2)
              end do !wi
            end do ! ii2

            !dldv Vg *********************************************************

            tr1=0
            do i=1,trn
              do j=1,trn
                tr1=tr1+pm(i,j)*pig(j,i)
              end do
            end do
            tr1=-0.5*tr1

            !tmpm=matmul(matmul(transpose(py),pig),py)
            call dgemm ('T','N',1,trn,trn,1.0D0,py,trn,pig,trn,0.0D0,py_tmp1,1)
            call dgemm ('N','N',1,1,trn,1.0D0,py_tmp1,1,py,trn,0.0D0,tmpm,1)

            dldv(idum+nbc,1)=tmpm(1,1)*0.5+tr1 !************************************
            !print*,'dldv:',tnk+nb*k-nb+xi,1,dldv(tn+nb*k-nb+xi,1),tn+nb*k-nb+xi
          end do !xi2
          end do !xi

          deallocate(phi,km2,tmp_wk,tmp_vcva,phi2,tmp_wk2)
        end do  ! jj
      end do  ! ii


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

  !write(41,*)'transformed ***'
  !do wi=1,mn
  !  do xi=1,tn
  !    write(41,'(a7,100f14.4)')'Va',vcva(wi,xi,xi)!,sqrt(sdmI(tn+wi*nb-nb+xi,tn+wi*nb-nb+xi))
  !    sum_v(xi)=sum_v(xi)+vcva(wi,xi,xi)
  !  end do
  !  nbc=0
  !  do xi=2,tn
  !    do yi=1,xi-1
  !      nbc=nbc+1
  !      write(41,'(a7,100f14.4)')'cov',vcva(wi,xi,yi)!,sqrt(sdmI(tn+wi*nb-nb+tn+nbc,tn+wi*nb-nb+tn+nbc))
  !    end do
  !  end do
  !  write(41,*)''
  !end do


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

if (fl11.ne.'null') then

  print*,'for BLUP solution after convergence *******************'

  !make it triangle -> full ****************
  do xi=1,tn
    do yi=1,xi
      do k=1,mn
        vcva(k,yi,xi)=vcva(k,xi,yi)
        !print*,vcva
      end do
    end do
  end do !***********************************

  do k=1,mn

    !allocate(gr_km_phi(sum(m_tnk(1:ttn,k)),tn))
    allocate(gr_km_phi(sum(m_tnk(1:ttn,k)),trn))

    do ii=1,ttn
      allocate(phi(rnm(ii),m_tnk(ii,k)),km(m_tnk(ii,k),m_tnk(ii,k)))
      allocate(tmp_wk(rnm(ii)),km_phi(m_tnk(ii,k),rnm(ii)))

      tmp_wk=m_wk(k,ii,1:rnm(ii))
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
        !allocate(phi(rnm(ii),m_tnk(ii,k)),km(m_tnk(ii,k),m_tnk(jj,k)))
        allocate(phi(rnm(jj),m_tnk(jj,k)),km(m_tnk(ii,k),m_tnk(jj,k)))
        !allocate(tmp_wk(rnm(ii)),km_phi(m_tnk(ii,k),rnm(jj)))
        allocate(tmp_wk(rnm(jj)),km_phi(m_tnk(ii,k),rnm(jj)))

        !tmp_wk=m_wk(ii,1:rnm(ii))
        tmp_wk=m_wk(k,jj,1:rnm(jj))
        !call pol (m_tnk(ii,k),tmp_wk,rnm(ii),phi)
        call pol (m_tnk(jj,k),tmp_wk,rnm(jj),phi)

        !allocate(phi2(rnm(jj),m_tnk(jj,k)),tmp_wk2(rnm(jj)))
        allocate(phi2(rnm(ii),m_tnk(ii,k)),tmp_wk2(rnm(ii)))
        !allocate(km_phi2(m_tnk(jj,k),rnm(ii)))
        allocate(km_phi2(m_tnk(jj,k),rnm(ii)))

        !tmp_wk2=m_wk(jj,1:rnm(jj))
        tmp_wk2=m_wk(k,ii,1:rnm(ii))
        !call pol (m_tnk(jj,k),tmp_wk2,rnm(jj),phi2)
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
          do j=1,rnm(yi)
            !v1=v1+mbin(k,yidx(yi,j),zi)*gr_km_phi(xi,yi)*py(trny-rnm(yi)+j,1)
            !v1=v1+mbin(k,yidx(yi,j),zi)*gr_km_phi(xi,j)*py(trny-rnm(yi)+j,1)
            v1=v1+mbin(k,yidx(yi,j),zi)*gr_km_phi(xi,trny-rnm(yi)+j)*py(trny-rnm(yi)+j,1)
          end do
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
    pv = huge ( pv )
    q = huge ( q )
    x = beta(xi,1)**2/xvmx(xi,xi)
    sd = 1.0D+00 !degree of greedom for cdfchi
    which=1  !Calculate P and Q from X, MEAN and SD;
    call cdfchi ( which, pv, q, x, sd, status, bound )
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

end subroutine

