
subroutine aireml_m_h2_dsgn (mbin,yidx,yv,rn,rnm,trn,pedn,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,fl24,xmm,nit,conv,wt_res,fl_dsgn,cov_idx)

!***********************************************************************
!aireml_m_h2_dsgn: multivariate additive + maternal genetic effects  
!S. Hong Lee (2020)

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
integer::yidx(tn,rn),rnm(tn),fnm(tn),cov_idx,meen,mee(rn)
integer::dsgn(rn,2),dsgn2(tn,rn,2)
!integer::rn1,rn2,yidx1(rn),yidx2(rn)
!double precision::mbin(mn,pedn,pedn),yv(rn,tn),tyv(trn,1)
real::mbin(mn,pedn,pedn),yv(rn,tn),sv_fact
double precision::resw(trn),resw2(pedn,tn),tyv(trn,1)

double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),conv
double precision::vcva(mn,tn,tn),vcve(tn,tn),vvc(tn)
double precision::fix_vcva(mn,tn,tn),fix_vcve(tn,tn)
double precision::vcva_rsv(mn,tn,tn),vcve_rsv(tn,tn)
double precision::blup_ebv(pedn,tn,mn),beta(tfn,1),blup_py(pedn,tn,mn) !,cor
double precision::blup_r(pedn,tn,mn),blup_mee(pedn,tn),blup_mee_r(pedn,tn)
double precision,allocatable::tmp_v(:,:)

integer::i,j,m,k,l,io
double precision::sum_h2,sum_v(tn)


double precision::LC,LA,LG,ypy,x1,x2,x3,x10,x11,v1,v2,v10,v11,v22,v222
double precision::y1,y2,y3,y10,y11,z1,z2,z3,z10,z11
double precision::LKHP,MLKH,mva(mn),mva2(mn),mve,mve2,mcov(mn)
double precision::h2(mn),h2n(mn),tr1,tr2(mn),tr3

double precision::xm(trn,tfn),xmm(tn,pedn,10000)

double precision::xvmx(tfn,tfn),xvm(tfn,trn) 
double precision::xvm2(tfn,trn),xvm3(trn,tfn) 
double precision::pm(trn,trn),py(trn,1),xb(trn,1)

double precision::v_tmp(trn),v2_tmp(trn)

double precision::aim(tn+(tn**2-tn)/2+mn*(tn+(tn**2-tn)/2),tn+(tn**2-tn)/2+mn*(tn+(tn**2-tn)/2))
double precision::sdmI(tn+(tn**2-tn)/2+mn*(tn+(tn**2-tn)/2),tn+(tn**2-tn)/2+mn*(tn+(tn**2-tn)/2))
double precision::up(tn+(tn**2-tn)/2+mn*(tn+(tn**2-tn)/2),1),dldv(tn+(tn**2-tn)/2+mn*(tn+(tn**2-tn)/2),1)



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
character(len=64)::fl4,fl24,fl5,fl11,fl11r,wt_res,cdum_t(tn),fl_dsgn
logical::file_exist

!CDF
real ( kind = 8 ) bound,mean,p,q,sd,x
integer ( kind = 4 ) status,which


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

!desing matrix reading and make Z2AZ2', (Z1AZ2' + Z2AZ1')
open (unit=52,file=fl_dsgn,status='old')
  do i=1,rn  !rn=pedn
    read(52,*)cdum,cdum,wi,wj
      dsgn(i,1)=wi
      dsgn(i,2)=wj
  end do
  close(52)

  do xi=1,tn
    k=0
    do i=1,rn
      if (yv(i,xi).ne.-99999) then
        k=k+1
        dsgn2(xi,k,1)=dsgn(i,1)
        dsgn2(xi,k,2)=dsgn(i,2)
        !print*,dsgn2(xi,k,:)
      end if
    end do
    if (rnm(xi).ne.k) then
      print*,k,rnm(xi),"check >>>"
    end if
  end do
  !pause

  mbin(2,:,:)=0
  do i=1,rn  !rn=pedn
  !do i=1,rnm(1)
    do j=1,i
      if (dsgn(i,2).ne.0 .and. dsgn(j,2).ne.0) then
        mbin(2,dsgn(i,1),dsgn(j,1))=mbin(1,dsgn(i,2),dsgn(j,2))
        mbin(2,dsgn(j,1),dsgn(i,1))=mbin(2,dsgn(i,1),dsgn(j,1))
        !print*,i,j,mbin(2,i,j)
        !pause
      end if
    end do
  end do
  !if (cov_idx==1) then ! with covariane between g and m too
  if (cov_idx >= 1) then ! with covariane between g and m too
    mbin(3,:,:)=0
    do i=1,rn  !rn=pedn
    !do i=1,rnm(1) 
      do j=1,i
        if (dsgn(i,1).ne.0 .and. dsgn(j,2).ne.0) then
          !mbin(3,i,j)=mbin(3,i,j)+mbin(1,dsgn(i,1),dsgn(j,2))
          mbin(3,dsgn(i,1),dsgn(j,1))=mbin(3,dsgn(i,1),dsgn(j,1))+mbin(1,dsgn(i,1),dsgn(j,2))
        end if
        if (dsgn(j,1).ne.0 .and. dsgn(i,2).ne.0) then
          !mbin(3,i,j)=mbin(3,i,j)+mbin(1,dsgn(j,1),dsgn(i,2))
          mbin(3,dsgn(i,1),dsgn(j,1))=mbin(3,dsgn(i,1),dsgn(j,1))+mbin(1,dsgn(j,1),dsgn(i,2))
        end if
        mbin(3,dsgn(j,1),dsgn(i,1))=mbin(3,dsgn(i,1),dsgn(j,1))
      end do
    end do

    if (cov_idx == 2) then ! with maternal environment
      mbin(4,:,:)=0
      do i=1,rn  !rn=pedn
      !do i=1,rnm(1)
        do j=1,i
          if (dsgn(i,2).ne.0 .and. dsgn(j,2).ne.0) then
            if (dsgn(i,2)==dsgn(j,2)) then
              mbin(4,dsgn(i,1),dsgn(j,1))=1
            end if
            mbin(4,dsgn(j,1),dsgn(i,1))=mbin(4,dsgn(i,1),dsgn(j,1))
          end if
        end do
      end do
    end if
  end if

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
!pause

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
        !print*,yv(i,xi),x1,x2
      end if
    end do
    !print*,x1,x2
    !pause
    vcve(xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))/2
    vcva(:,xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))/(2*mn)
  end do
  !print*,vcve(1,1)
  
  inquire(file=fl4,exist=file_exist)
  if (file_exist) then
    open (UNIT=44,FILE=fl4,STATUS='old')
    do i=1,tn
      read(44,*,iostat=io)cdum,vcve(i,i)
      if (io.ne.0) exit
    end do
    do xi=1,tn
      do yi=1,xi-1
        read(44,*,iostat=io)cdum,vcve(xi,yi)
        if (io.ne.0) exit
      end do
    end do
    !print*,vcve

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

  inquire(file=fl24,exist=file_exist)
  if (file_exist) then
    open (UNIT=51,FILE=fl24,STATUS='old')
    do i=1,tn
      read(51,*,iostat=io)cdum,fix_vcve(i,i)
      if (io.ne.0) exit
    end do
    do xi=1,tn
      do yi=1,xi-1
        read(51,*,iostat=io)cdum,fix_vcve(xi,yi)
        if (io.ne.0) exit
      end do
    end do
    !print*,vcve

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
    do xi=1,tn
      do yi=1,xi-1
        if (fix_vcve(xi,yi) .ne. -99999) vcve(xi,yi)=fix_vcve(xi,yi)
      end do
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
  !*****************************
  

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
          if (yidx(xi,i)==yidx(yi,j)) then
            pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)+vcve(xi,yi)*sqrt(resw(trnx-rnm(xi)+i)*resw(trny-rnm(yi)+j)) 
          end if
          if (xi.ne.yi) then
            pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)
          end if
        end do
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

    !BLUE = (X'VI X)I XVI y **********************
    call dgemm ('T','N',trn,tfn,tfn,1.0D0,xvm,tfn,xvmx,tfn,0.0D0,xvm3,trn)
    call dgemm ('T','N',tfn,1,trn,1.0D0,xvm3,trn,tyv,trn,0.0D0,beta,tfn)

    call dgemm('N','N',trn, trn, tfn, -1.0D0, xvm3, trn, xvm, tfn, 1.0D0,pm, trn)

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
    !pause

    call cpu_time (t2_cpu)
    !print*,'LKH:',t2_cpu-t1_cpu

    !if ((LKH>-99999999 .and. LKH<99999999) .and. (LKH>=LKHP)) then
    if ((LKH>-99999999 .and. LKH<99999999) ) then
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


    PRINT '(a7,100f14.4)','LKH',LKH
    do xi=1,tn
      PRINT '(a7,100f14.4)','Ve',vcve(xi,xi)
    end do
    do xi=2,tn
      do yi=1,xi-1
        PRINT '(a7,100f14.4)','cov',vcve(xi,yi)
      end do
    end do
    print*,''
    !do wi=1,mn
    do wi=1,2
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
    if (mn>=3) then
      do xi=1,tn
        PRINT '(a7,100f14.4)','covgm',vcva(3,xi,xi)
      end do
      do xi=2,tn
        do yi=1,xi-1
          PRINT '(a7,100f14.4)','cov',vcva(3,xi,yi)
        end do
      end do
      print*,''
      if (mn==4) then
        do xi=1,tn
          PRINT '(a7,100f14.4)','Va',vcva(4,xi,xi)
        end do
        do xi=2,tn
          do yi=1,xi-1
            PRINT '(a7,100f14.4)','cov',vcva(4,xi,yi)
          end do
        end do
        print*,''
      end if
      !pause
    end if


    ! AI matrix
    py=MATMUL(pm,tyv)       !for py

    !print*,'ai matrix'
    call cpu_time (t1_cpu)
    aim=0
    dldv=0

    nb=tn+(tn**2-tn)/2  !no block for random effects (Vg, cov)
    !print*,nb,tn
    !Ve (diagonal)  **********************************
    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)   !total rn upto xi
      v2=0
      ! $OMP parallel do simd reduction(+:v2)
      do i=1,rnm(xi)
        v1=0
        do j=1,rnm(xi)
          v1=v1+py(trnx-rnm(xi)+j,1)*resw(trnx-rnm(xi)+j)*resw(trnx-rnm(xi)+i)*pm(trnx-rnm(xi)+j,trnx-rnm(xi)+i)
        end do
        v2=v2+v1*py(trnx-rnm(xi)+i,1)
      end do
      ! $OMP end parallel do
      aim(xi,xi)=v2
    end do

    !dldv Ve ***********************************************************
    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)   !total rn upto xi
      tr1=0
      ! $OMP parallel do simd reduction(+:tr1)
      do i=1,rnm(xi)
        tr1=tr1+resw(trnx-rnm(xi)+i)*pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)
      end do
      ! $OMP end parallel do simd
      tr1=-0.5*tr1
      v2=0
      ! $OMP parallel do simd reduction(+:v2)
      do i=1,rnm(xi)
        v2=v2+py(trnx-rnm(xi)+i,1)*resw(trnx-rnm(xi)+i)*py(trnx-rnm(xi)+i,1)
      end do
      ! $OMP end parallel do simd 
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
        ! $OMP parallel do simd reduction(+:v2)
        do i=1,rnm(xi)
          v1=0
          ! $OMP parallel do simd reduction(+:v1)
          do j=1,rnm(yi)
            !v1=v1+py(trny-rnm(yi)+j,1)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
            v1=v1+py(trny-rnm(yi)+j,1)*resw(trny-rnm(yi)+j)*resw(trnx-rnm(xi)+i)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
          end do
          ! $OMP end parallel do simd
          v2=v2+v1*py(trnx-rnm(xi)+i,1)
        end do 
        ! $OMP end parallel do simd
        aim(xi,yi)=v2
        !print*,v2
      end do
    end do
    !pause

    !for cove (diagonal) *************************************************
    nbc=0
    trnx=rnm(1)
    do xi=2,tn
      trnx=trnx+rnm(xi)
      trny=0
      do yi=1,xi-1
        v_tmp=0
        nbc=nbc+1
        trny=trny+rnm(yi)
        v2=0
        ! $OMP parallel do simd reduction(+:v2)
        do i=1,rnm(xi)
          v10=0
          ! $OMP parallel do simd reduction(+:v10)
          do j=1,rnm(yi)
            if (yidx(xi,i)==yidx(yi,j)) then
              v10=v10+py(trny-rnm(yi)+j,1)*resw(trny-rnm(yi)+j) !py*pig
            end if
          end do
          ! $OMP end parallel do simd
          v_tmp(trnx-rnm(xi)+i)=v10
        end do
        ! $OMP end parallel do simd
        !$OMP PARALLEL DO PRIVATE (i,j,v10)
        do i=1,rnm(yi)
          v10=0
          do j=1,rnm(xi)
            if (yidx(xi,j)==yidx(yi,i)) then
              v10=v10+py(trnx-rnm(xi)+j,1)*resw(trnx-rnm(xi)+j)
            end if
          end do
          v_tmp(trny-rnm(yi)+i)=v10
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO PRIVATE (i,j,v10)
        do i=1,trn     !check (may not trn)
          v10=0
          do j=1,rnm(xi)
            v10=v10+v_tmp(trnx-rnm(xi)+j)*pm(trnx-rnm(xi)+j,i) !*pm
          end do
          do j=1,rnm(yi)
            v10=v10+v_tmp(trny-rnm(yi)+j)*pm(trny-rnm(yi)+j,i) !*pm
          end do
          v2_tmp(i)=v10
        end do
        !$OMP END PARALLEL DO

        ! $OMP parallel do 
        do i=1,rnm(xi)
          v1=0
          do j=1,rnm(yi)
            if (yidx(xi,i)==yidx(yi,j)) then
              v1=v1+v2_tmp(trny-rnm(yi)+j)*resw(trny-rnm(yi)+j) !*pig
            end if
          end do
          v2=v2+v1*py(trnx-rnm(xi)+i,1)                                !*py
        end do
        ! $OMP end parallel do 

        ! $OMP parallel do 
        do i=1,rnm(yi)
          v1=0
          do j=1,rnm(xi)
            if (yidx(xi,j)==yidx(yi,i)) then
              v1=v1+v2_tmp(trnx-rnm(xi)+j)*resw(trnx-rnm(xi)+j) !*pig
            end if
          end do
          v2=v2+v1*py(trny-rnm(yi)+i,1)                                !*py
        end do
        ! $OMP end parallel do 
        aim(tn+nbc,tn+nbc)=v2     !check

        !for cove x ve (off diagonal) *****************************************
        trnv=0
        do vi=1,tn
          trnv=trnv+rnm(vi)
          v2=0
          do i=1,rnm(vi)
            !v2=v2+v2_tmp(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)   !*py
            v2=v2+v2_tmp(trnv-rnm(vi)+i)*resw(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)   !*py
          end do
          aim(tn+nbc,vi)=v2          !check with Ve1
        end do

        !for cove x cove (off diagonal) !**************************************
        trnv=rnm(1)
        nbc2=0    !no block for cov (vi,ui)
        do vi=2,tn
          trnv=trnv+rnm(vi)
          trnu=0
          ! $OMP parallel do simd reduction(+:trnu)
          do ui=1,(vi-1)
            nbc2=nbc2+1

            if (nbc2<nbc) then  !off diagonal within cove ***************
              trnu=trnu+rnm(ui)
              v2=0
              do i=1,rnm(vi)
                v1=0
                do j=1,rnm(ui)
                  if (yidx(ui,j)==yidx(vi,i)) then
                    !v1=v1+v2_tmp(trnu-rnm(ui)+j) !*pig_tmp
                    v1=v1+v2_tmp(trnu-rnm(ui)+j)*sqrt(resw(trnu-rnm(ui)+j)*resw(trnv-rnm(vi)+i)) !*pig_tmp
                  end if
                end do
                v2=v2+v1*py(trnv-rnm(vi)+i,1) !*py
              end do
              do i=1,rnm(ui)
                v1=0
                do j=1,rnm(vi)
                  if (yidx(ui,i)==yidx(vi,j)) then
                    !v1=v1+v2_tmp(trnv-rnm(vi)+j) !*pig_tmp
                    v1=v1+v2_tmp(trnv-rnm(vi)+j)*sqrt(resw(trnv-rnm(vi)+j)*resw(trnu-rnm(ui)+i)) !*pig_tmp
                  end if
                end do
                v2=v2+v1*py(trnu-rnm(ui)+i,1)      
              end do
              aim(tn+nbc,tn+nbc2)=v2
              !print*,v2
            end if ! *****************************************************
          end do !ui
          ! $OMP end parallel do simd
        end do !vi

      end do  !yi
    end do  !xi

    !dldv for cove ******************************************************
    nbc=0     !no block for cov (xi, yi)
    trnx=rnm(1)
    do xi=2,tn
      trnx=trnx+rnm(xi)
      trny=0
      ! $OMP parallel do simd reduction(+:trny)
      do yi=1,xi-1
        nbc=nbc+1
        trny=trny+rnm(yi)
        tr1=0
        do i=1,rnm(xi)
          v10=0
          do j=1,rnm(yi)
            if (yidx(yi,j)==yidx(xi,i)) then
              !v10=v10+pm(trnx-rnm(xi)+i,trny-rnm(yi)+j) !pm*pig
              v10=v10+pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)*sqrt(resw(trnx-rnm(xi)+i)*resw(trny-rnm(yi)+j)) !pm*pig
            end if
          end do
          tr1=tr1+v10
        end do
        do i=1,rnm(yi)
          v10=0
          do j=1,rnm(xi)
            if (yidx(xi,j)==yidx(yi,i)) then
              !v10=v10+pm(trny-rnm(yi)+i,trnx-rnm(xi)+j) !pm*pig
              v10=v10+pm(trny-rnm(yi)+i,trnx-rnm(xi)+j)*sqrt(resw(trny-rnm(yi)+i)*resw(trnx-rnm(xi)+j)) !pm*pig
            end if
          end do
          tr1=tr1+v10
        end do
        tr1=-0.5*tr1

        v2=0
        do i=1,rnm(xi)
          v10=0
          do j=1,rnm(yi)
            if (yidx(xi,i)==yidx(yi,j)) then
              !v10=v10+py(trny-rnm(yi)+j,1) !py*pig
              v10=v10+py(trny-rnm(yi)+j,1)*resw(trny-rnm(yi)+j) !py*pig
            end if
          end do
          v2=v2+v10*py(trnx-rnm(xi)+i,1)
        end do
        do i=1,rnm(yi)
          v10=0
          do j=1,rnm(xi)
            if (yidx(xi,j)==yidx(yi,i)) then
              !v10=v10+py(trnx-rnm(xi)+j,1) !py*pig
              v10=v10+py(trnx-rnm(xi)+j,1)*resw(trnx-rnm(xi)+j) !py*pig
            end if
          end do
          v2=v2+v10*py(trny-rnm(yi)+i,1)
        end do
        dldv(tn+nbc,1)=v2*0.5+tr1
        !print*,v2,tr1
      end do !yi
      ! $OMP end parallel do simd
    end do !xi



    !print*,'Vg'
    !for Vg
    do k=1,mn !no. random effects

      !for Vg (doagonal) ******************************************************* 
      trnx=0
      do xi=1,tn
        trnx=trnx+rnm(xi)
        v2=0
        ! $OMP parallel do simd reduction(+:v2)
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
        ! $OMP end parallel do simd

        aim(nb+nb*k-nb+xi,nb+nb*k-nb+xi)=v2

        !for Vg x Ve (off diagonal) !******************************************
        trnv=0
        do vi=1,tn
          trnv=trnv+rnm(vi)
          v2=0
          do i=1,rnm(vi)
            v2=v2+v2_tmp(trnv-rnm(vi)+i)*resw(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)
!*py
          end do
          aim(nb+nb*k-nb+xi,vi)=v2
        end do

        !for Vg x cove ***************************************************
        nbc=0  !for no. block for cove, i.e. 1:(2,1), 2:(3,1), 3:(3,2), ... 
        trnv=rnm(1)
        do vi=2,tn
          trnv=trnv+rnm(vi)   !total rn upto vi
          trnu=0
          do ui=1,(vi-1)
            nbc=nbc+1
            trnu=trnu+rnm(ui)
            v2=0
            ! $OMP parallel do simd reduction(+:v2)
            do i=1,rnm(vi)
              v1=0
              do j=1,rnm(ui)
                if (yidx(ui,j)==yidx(vi,i)) then
                  !v1=v1+v2_tmp(trnu-rnm(ui)+j) !*pig_tmp
                  v1=v1+v2_tmp(trnu-rnm(ui)+j)*sqrt(resw(trnu-rnm(ui)+j)*resw(trnv-rnm(vi)+i)) !*pig_tmp
                end if
              end do
              v2=v2+v1*py(trnv-rnm(vi)+i,1) !*py
            end do
            do i=1,rnm(ui)
              v1=0
              do j=1,rnm(vi)
                if (yidx(ui,i)==yidx(vi,j)) then
                  !v1=v1+v2_tmp(trnv-rnm(vi)+j) !*pig_tmp
                  v1=v1+v2_tmp(trnv-rnm(vi)+j)*sqrt(resw(trnv-rnm(vi)+j)*resw(trnu-rnm(ui)+i)) !*pig_tmp
                end if
              end do
              v2=v2+v1*py(trnu-rnm(ui)+i,1) !*py
            end do
            ! $OMP end parallel do simd
            aim(nb+nb*k-nb+xi,tn+nbc)=v2  !check
          end do !vi
        end do !ui


        !for Vg(1~[i-1]) x Vg(i) **********************************************
        trnv=0
        do vi=1,xi-1
          trnv=trnv+rnm(vi)
          v2=0
          ! $OMP parallel do simd reduction(+:v2)
          do i=1,rnm(vi)
            v1=0
            do j=1,rnm(vi)
              v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(k,yidx(vi,j),yidx(vi,i))   !*pig_tmp
            end do
            v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
          end do
          ! $OMP end parallel do simd
          aim(nb+nb*k-nb+xi,nb+nb*k-nb+vi)=v2
        end do


        !off diagonal *********************************************************
        do wi=1,k-1
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            ! $OMP parallel do simd reduction(+:v2)
            do i=1,rnm(vi)
              v1=0
              do j=1,rnm(vi)
                v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(wi,yidx(vi,j),yidx(vi,i))     !*pig_tmp
              end do
              v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
            end do
            ! $OMP end parallel do simd
            aim(nb+nb*k-nb+xi,nb+nb*wi-nb+vi)=v2   !with the xi th Vg for first trait
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
              ! $OMP parallel do simd reduction(+:v2)
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
              ! $OMP end parallel do simd
              !aim(2+k*3-3+xi,2+wi*3)=v2  !check
              aim(nb+nb*k-nb+xi,nb+nb*wi-nb+tn+nbc)=v2  !check
            end do

          end do

        end do !wi

      end do  !xi
      !pause


      !dldv Vg *********************************************************
      trnx=0
      ! $OMP parallel do simd reduction(+:trnx)
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
        dldv(nb+nb*k-nb+xi,1)=v2*0.5+tr1 !************************************
      end do
      ! $OMP end parallel do simd


      !for cova (diagonal) *************************************************
      nbc=0  !for no. block for covariance, i.e. 1: (2,1), 2: (3,1), 3:(3,2), ... 
      trnx=rnm(1)
      do xi=2,tn
        trnx=trnx+rnm(xi)
        trny=0
        ! $OMP parallel do simd reduction(+:trny)
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
          aim(nb+nb*k-nb+tn+nbc,nb+nb*k-nb+tn+nbc)=v2     !check

          !for cov x Ve (off diagonal) ************************************
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              !v2=v2+v2_tmp(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)         !*py
              v2=v2+v2_tmp(trnv-rnm(vi)+i)*resw(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)
!*py
            end do 
            aim(nb+nb*k-nb+tn+nbc,vi)=v2          !check with Ve1
          end do

          !cova x cove **********************************************
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
                  if (yidx(ui,j)==yidx(vi,i)) then
                    !v1=v1+v2_tmp(trnu-rnm(ui)+j) !*pig_tmp
                    v1=v1+v2_tmp(trnu-rnm(ui)+j)*sqrt(resw(trnu-rnm(ui)+j)*resw(trnv-rnm(vi)+i)) !*pig_tmp
                  end if
                end do
                v2=v2+v1*py(trnv-rnm(vi)+i,1) !*py
              end do
              do i=1,rnm(ui)
                v1=0
                do j=1,rnm(vi)
                  if (yidx(ui,i)==yidx(vi,j)) then
                    !v1=v1+v2_tmp(trnv-rnm(vi)+j) !*pig_tmp
                    v1=v1+v2_tmp(trnv-rnm(vi)+j)*sqrt(resw(trnv-rnm(vi)+j)*resw(trnu-rnm(ui)+i)) !*pig_tmp
                  end if
                end do
                v2=v2+v1*py(trnu-rnm(ui)+i,1)  !*py
              end do
              aim(nb+nb*k-nb+tn+nbc,tn+nbc2)=v2 
            end do !ui
          end do !vi


          !for cova x Vg **************************************************
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
            aim(nb+nb*k-nb+tn+nbc,nb+nb*k-nb+vi)=v2  !with the same Vg1 for first trait
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
                  v1=v1+v2_tmp(trnu-rnm(ui)+j)*mbin(k,yidx(ui,j),yidx(vi,i)) !*pig_tmp
                end do
                v2=v2+v1*py(trnv-rnm(vi)+i,1) !*py
              end do
              do i=1,rnm(ui)
                v1=0
                do j=1,rnm(vi)
                  v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(k,yidx(vi,j),yidx(ui,i)) !*pig_tmp
                end do
                v2=v2+v1*py(trnu-rnm(ui)+i,1) !*py
              end do
              aim(nb+nb*k-nb+tn+nbc,nb+nb*k-nb+tn+nbc2)=v2  

              end if ! ***********************************************************
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
              aim(nb+nb*k-nb+tn+nbc,nb+nb*wi-nb+vi)=v2  !with the wi th Vg for first trait
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
                aim(nb+nb*k-nb+tn+nbc,nb+nb*wi-nb+tn+nbc2)=v2  !with the wi th cov for first trait
              end do !ui
            end do !vi
          end do !wi
        end do !yi
        ! $OMP end parallel do simd

      end do !xi **********************************************

      !dldv ***************************************************************
      nbc=0     !no block for cov (xi, yi)
      trnx=rnm(1)
      do xi=2,tn
        trnx=trnx+rnm(xi)
        trny=0
        ! $OMP parallel do simd reduction(+:trny)
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
          dldv(nb+nb*k-nb+tn+nbc,1)=v2*0.5+tr1
        end do !yi
        ! $OMP end parallel do simd
      end do !xi


    end do !k

    call cpu_time (t2_cpu)
    print*,'derivatives done - time:',real(t2_cpu-t1_cpu)
    print*,''
 
    !do i=1,(2+mn*3)
    do i=1,nb+nb*mn
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
    do xi=1,tn
      do yi=1,xi-1
        wj=wj+1
        !print*,wj
        if (fix_vcve(xi,yi) .ne. -99999) then
            aim(:,wj)=0
            aim(wj,:)=0
          aim(wj,wj)=1
        end if
      end do
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
  

    !do i=1,nb+nb*mn
    !  print*,real(aim(i,:))
    !end do
    !print*,''
    !pause
    call cholesky_inv (aim,(nb+nb*mn),x1)

    !print*,x1
    !do i=1,nb+nb*mn
    !  print*,real(aim(i,:))
    !end do
    !print*,''
    !pause

    if (LKH-MLKH.gt.0) then
      MLKH=LKH
      mva=va
      mva2=va2
      mve=ve
      mve2=ve2
      mcov=cov
      sdmI=aim
    end if

    !if (ABS(LKH-LKHP).lt.0.001) goto 1000
    if (ABS(LKH-LKHP).lt.conv) goto 1000

    up=MATMUL(aim,dldv)
    !do i=1,nb+nb*mn
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
    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        vcve(xi,yi)=vcve(xi,yi)+up(tn+nbc,1)
      end do
    end do

    do wi=1,mn
      do xi=1,tn
        !print*,nb,wi,xi,nb+nb*wi-nb+xi
        vcva(wi,xi,xi)=vcva(wi,xi,xi)+up(nb+nb*wi-nb+xi,1)
        !constrain out of parameter space
        !if (vcva(wi,xi,xi).lt.0) vcva(wi,xi,xi)=0.001
      end do
      nbc=0
      do xi=2,tn
        do yi=1,xi-1
          nbc=nbc+1
          vcva(wi,xi,yi)=vcva(wi,xi,yi)+up(nb+nb*wi-nb+tn+nbc,1)
          !constrain out of parameter space
          !if (vcva(wi,xi,yi)/sqrt(vcva(wi,xi,xi)*vcva(wi,yi,yi))>=1) then
          !  vcva(wi,xi,yi)=sqrt(vcva(wi,xi,xi)*vcva(wi,yi,yi))
          !end if


        end do
      end do
    end do

    !print*,ve,ve2,va,va2,cov
    !print*,vcve,vcva

  END do !zi



  !PRINT*,'LKH not converged na na na'
  !WRITE(41,*)'LKH not converged na na na'
  va=mva;va2=mva2;ve=mve;ve2=mve2;LKH=MLKH;cov=mcov


1000 continue 


  open (UNIT=41,FILE=fl5,STATUS='unknown')

    !fixing some parameters
    wj=0
    do i=1,tn
      wj=wj+1
      if (fix_vcve(i,i) .ne. -99999) then
        sdmI(wj,wj)=0
      end if
    end do
    do xi=1,tn
      do yi=1,xi-1
        wj=wj+1
        if (fix_vcve(xi,yi) .ne. -99999) then
          sdmI(wj,wj)=0
        end if
      end do
    end do
    do wi=1,mn
      do xi=1,tn
        wj=wj+1
        if (fix_vcva(wi,xi,xi).ne.-99999) then
          sdmI(wj,wj)=0
        end if
      end do
      do xi=1,tn
        do yi=1,xi-1
          wj=wj+1
          if (fix_vcva(wi,xi,yi).ne.-99999) then
            sdmI(wj,wj)=0
          end if
        end do
      end do
    end do
    !*****************************

  do xi=1,tn
    write(41,'(a7,100f14.4)')'Ve',vcve(xi,xi),sqrt(sdmI(xi,xi))
  end do
  nbc=0
  do xi=2,tn
    do yi=1,xi-1
      nbc=nbc+1
      write(41,'(a7,100f14.4)')'cove',vcve(xi,yi),sqrt(sdmI(tn+nbc,tn+nbc))
    end do
  end do
  write(41,*)''

  do xi=1,tn
    sum_v(xi)=vcve(xi,xi)
  end do
  !do wi=1,mn
  do wi=1,2
    do xi=1,tn
      write(41,'(a7,100f14.4)')'Va',vcva(wi,xi,xi),sqrt(sdmI(nb+wi*nb-nb+xi,nb+wi*nb-nb+xi))
      sum_v(xi)=sum_v(xi)+vcva(wi,xi,xi)
    end do
    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        write(41,'(a7,100f14.4)')'cov',vcva(wi,xi,yi),sqrt(sdmI(nb+wi*nb-nb+tn+nbc,nb+wi*nb-nb+tn+nbc))
      end do
    end do
    write(41,*)''
  end do
  if (mn>=3) then
    do xi=1,tn
      write(41,'(a7,100f14.4)')'covgm',vcva(3,xi,xi),sqrt(sdmI(nb+3*nb-nb+xi,nb+3*nb-nb+xi))
      sum_v(xi)=sum_v(xi)+vcva(3,xi,xi)*2 !single trait >>> check for multi-trait
    end do
    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        write(41,'(a7,100f14.4)')'cov',vcva(3,xi,yi),sqrt(sdmI(nb+3*nb-nb+tn+nbc,nb+3*nb-nb+tn+nbc))
      end do
    end do
    write(41,*)''

    if (mn==4) then
      do xi=1,tn
        write(41,'(a7,100f14.4)')'Va',vcva(4,xi,xi),sqrt(sdmI(nb+4*nb-nb+xi,nb+4*nb-nb+xi))
        sum_v(xi)=sum_v(xi)+vcva(4,xi,xi)  !single trait >>> check for multi-trait
      end do
      nbc=0
      do xi=2,tn
        do yi=1,xi-1
          nbc=nbc+1
          write(41,'(a7,100f14.4)')'cov',vcva(4,xi,yi),sqrt(sdmI(nb+4*nb-nb+tn+nbc,nb+4*nb-nb+tn+nbc))
        end do
      end do
      write(41,*)''
    end if
    write(41,*)''
  end if

  !total Var(V) *********************
  do xi=1,tn
    vvc(xi)=sdmI(xi,xi)
    !do wi=1,mn
    do wi=1,2
      vvc(xi)=vvc(xi)+sdmI(nb+nb*wi-nb+xi,xi)*2
      !do wj=1,mn
      do wj=1,2
        vvc(xi)=vvc(xi)+sdmI(nb+nb*wi-nb+xi,nb+nb*wj-nb+xi)
      end do
      if (mn>=3) then
        vvc(xi)=vvc(xi)+sdmI(nb+nb*wi-nb+xi,nb+nb*3-nb+xi)*2
        if (mn==4) then
          vvc(xi)=vvc(xi)+sdmI(nb+nb*wi-nb+xi,nb+nb*4-nb+xi)
        end if
      end if
    end do
    if (mn >= 3) then
      vvc(xi)=vvc(xi)+sdmI(nb+nb*3-nb+xi,xi)*4
      do wj=1,2
        vvc(xi)=vvc(xi)+sdmI(nb+nb*3-nb+xi,nb+nb*wj-nb+xi)*2
      end do
      vvc(xi)=vvc(xi)+sdmI(nb+nb*3-nb+xi,nb+nb*3-nb+xi)*2
      if (mn == 4) then
        vvc(xi)=vvc(xi)+sdmI(nb+nb*4-nb+xi,xi)*2
        do wj=1,2
          vvc(xi)=vvc(xi)+sdmI(nb+nb*4-nb+xi,nb+nb*wj-nb+xi)
        end do
        vvc(xi)=vvc(xi)+sdmI(nb+nb*4-nb+xi,nb+nb*3-nb+xi)*2
        vvc(xi)=vvc(xi)+sdmI(nb+nb*4-nb+xi,nb+nb*4-nb+xi)
      end if
    end if

  end do

  !print*,vvc

  !do wi=1,mn
  do wi=1,2
    do xi=1,tn
      x1=vcva(wi,xi,xi)/sum_v(xi)
      x11=sdmI(nb+nb*wi-nb+xi,xi) !for the row ************

      !do wj=1,mn
      do wj=1,2
        x11=x11+sdmI(nb+nb*wi-nb+xi,nb+nb*wj-nb+xi)
      end do
      if (mn>=3) then
        x11=x11+sdmI(nb+nb*wi-nb+xi,nb+nb*3-nb+xi)*2
        if (mn==4) then
          x11=x11+sdmI(nb+nb*wi-nb+xi,nb+nb*4-nb+xi)
        end if
      end if
      !print*,x11

      !SE **********************************
      !x2=x1**2*(sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+xi)/vcva(wi,xi,xi)**2+vvc(xi)/sum_v(xi)**2-2*x11/(vcva(wi,xi,xi)*sum_v(xi)))
      x2=x1**2*(sdmI(nb+nb*wi-nb+xi,nb+nb*wi-nb+xi)/vcva(wi,xi,xi)**2+vvc(xi)/sum_v(xi)**2-2*x11/(vcva(wi,xi,xi)*sum_v(xi)))
      write(41,'(a7,100f14.4)')'h2',x1,sqrt(x2)

    end do
    write(41,*)''

    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        z1=vcva(wi,xi,yi)/sqrt(vcva(wi,xi,xi)*vcva(wi,yi,yi))
        !var(cor) *************************************
        z2=sdmI(nb+nb*wi-nb+xi,nb+nb*wi-nb+xi)/(4*vcva(wi,xi,xi)**2)
        z2=z2+sdmI(nb+nb*wi-nb+yi,nb+nb*wi-nb+yi)/(4*vcva(wi,yi,yi)**2)
        z2=z2+sdmI(nb+nb*wi-nb+tn+nbc,nb+nb*wi-nb+tn+nbc)/(vcva(wi,xi,yi)**2)
        z2=z2+2*sdmI(nb+nb*wi-nb+xi,nb+nb*wi-nb+yi)/(4*vcva(wi,xi,xi)*vcva(wi,yi,yi))
        z2=z2-2*sdmI(nb+nb*wi-nb+xi,nb+nb*wi-nb+tn+nbc)/(2*vcva(wi,xi,xi)*vcva(wi,xi,yi))
        z2=z2-2*sdmI(nb+nb*wi-nb+tn+nbc,nb+nb*wi-nb+yi)/(2*vcva(wi,xi,yi)*vcva(wi,yi,yi))
        z2=z2*(z1**2)
        write(41,'(a7,100f14.4)')'cor',z1,sqrt(z2)
      end do
    end do
    write(41,*)''
  end do
  if (mn>=3) then
    do xi=1,tn
      z1=vcva(3,xi,xi)/sqrt(vcva(1,xi,xi)*vcva(2,xi,xi))
      !var(cor) *************************************

      z2=sdmI(nb+nb*1-nb+xi,nb+nb*1-nb+xi)/(4*vcva(1,xi,xi)**2)
      z2=z2+sdmI(nb+nb*2-nb+xi,nb+nb*2-nb+xi)/(4*vcva(2,xi,xi)**2)
      z2=z2+sdmI(nb+nb*3-nb+xi,nb+nb*3-nb+xi)/(vcva(3,1,1)**2)
      z2=z2+2*sdmI(nb+nb*1-nb+xi,nb+nb*2-nb+xi)/(4*vcva(1,xi,xi)*vcva(2,xi,xi))
      z2=z2-2*sdmI(nb+nb*1-nb+xi,nb+nb*3-nb+xi)/(2*vcva(1,xi,xi)*vcva(3,xi,xi))
      z2=z2-2*sdmI(nb+nb*3-nb+xi,nb+nb*2-nb+xi)/(2*vcva(3,xi,xi)*vcva(2,xi,xi))
      z2=z2*(z1**2)
      write(41,'(a7,100f14.4)')'corgm',z1,sqrt(z2)
    end do
    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        z1=vcva(3,xi,yi)/sqrt(vcva(3,xi,xi)*vcva(3,yi,yi))
        !var(cor) *************************************
        z2=sdmI(nb+nb*3-nb+xi,nb+nb*3-nb+xi)/(4*vcva(3,xi,xi)**2)
        z2=z2+sdmI(nb+nb*3-nb+yi,nb+nb*3-nb+yi)/(4*vcva(3,yi,yi)**2)
        z2=z2+sdmI(nb+nb*3-nb+tn+nbc,nb+nb*3-nb+tn+nbc)/(vcva(3,xi,yi)**2)
        z2=z2+2*sdmI(nb+nb*3-nb+xi,nb+nb*3-nb+yi)/(4*vcva(3,xi,xi)*vcva(3,yi,yi))
        z2=z2-2*sdmI(nb+nb*3-nb+xi,nb+nb*3-nb+tn+nbc)/(2*vcva(3,xi,xi)*vcva(3,xi,yi))
        z2=z2-2*sdmI(nb+nb*3-nb+tn+nbc,nb+nb*3-nb+yi)/(2*vcva(3,xi,yi)*vcva(3,yi,yi))
        z2=z2*(z1**2)
        write(41,'(a7,100f14.4)')'cor',z1,sqrt(z2)
      end do
    end do
    write(41,*)''

    if (mn==4) then
      do xi=1,tn
        x1=vcva(4,xi,xi)/sum_v(xi)
        x11=sdmI(nb+nb*4-nb+xi,xi) !for the row ************
        x11=x11+sdmI(nb+nb*4-nb+xi,nb+nb*1-nb+xi)
        x11=x11+sdmI(nb+nb*4-nb+xi,nb+nb*2-nb+xi)
        x11=x11+sdmI(nb+nb*4-nb+xi,nb+nb*3-nb+xi)*2
        x11=x11+sdmI(nb+nb*4-nb+xi,nb+nb*4-nb+xi)
        x2=x1**2*(sdmI(nb+nb*4-nb+xi,nb+nb*4-nb+xi)/vcva(4,xi,xi)**2+vvc(xi)/sum_v(xi)**2-2*x11/(vcva(4,xi,xi)*sum_v(xi)))
        write(41,'(a7,100f14.4)')'h2',x1,sqrt(x2)
      end do
      nbc=0
      do xi=2,tn
        do yi=1,xi-1
          nbc=nbc+1
          z1=vcva(4,xi,yi)/sqrt(vcva(4,xi,xi)*vcva(4,yi,yi))
          !var(cor) *************************************
          z2=sdmI(nb+nb*4-nb+xi,nb+nb*4-nb+xi)/(4*vcva(4,xi,xi)**2)
          z2=z2+sdmI(nb+nb*4-nb+yi,nb+nb*4-nb+yi)/(4*vcva(4,yi,yi)**2)
          z2=z2+sdmI(nb+nb*4-nb+tn+nbc,nb+nb*4-nb+tn+nbc)/(vcva(4,xi,yi)**2)
          z2=z2+2*sdmI(nb+nb*4-nb+xi,nb+nb*4-nb+yi)/(4*vcva(4,xi,xi)*vcva(4,yi,yi))
          z2=z2-2*sdmI(nb+nb*4-nb+xi,nb+nb*4-nb+tn+nbc)/(2*vcva(4,xi,xi)*vcva(4,xi,yi))
          z2=z2-2*sdmI(nb+nb*4-nb+tn+nbc,nb+nb*4-nb+yi)/(2*vcva(4,xi,yi)*vcva(4,yi,yi))
          z2=z2*(z1**2)
          write(41,'(a7,100f14.4)')'cor',z1,sqrt(z2)
        end do
      end do
      write(41,*)''
    end if
  end if
  write(41,*)''

  write(41,*)''
  write(41,'(a7,100f14.4)')'LKH',LKH

  write(41,*)''
  do i=1,(nb+nb*mn)
    !write(41,*)real(sdmI(i,1:i))
    !write(41,'(100f12.5)')real(sdmI(i,1:i))
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

    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)
      do zi=1,pedn
        !do k=1,mn
        do k=1,2
          v1=0;v2=0
          trny=0
          do yi=1,tn
            trny=trny+rnm(yi)
            do j=1,rnm(yi)
              if (dsgn2(yi,j,k).ne.0) then
                v1=v1+mbin(1,dsgn2(yi,j,k),zi)*vcva(k,xi,yi)*py(trny-rnm(yi)+j,1)
              end if
              if (dsgn2(yi,j,(3-k)).ne.0) then
                v1=v1+mbin(1,dsgn2(yi,j,(3-k)),zi)*vcva(3,xi,yi)*py(trny-rnm(yi)+j,1)
              end if
              if (dsgn2(yi,j,k)==zi) then
                v2=v2+vcva(k,xi,yi)*py(trny-rnm(yi)+j,1)  !not sure for this >>> check
              end if
            end do
          end do
          blup_ebv(zi,xi,k)=v1
          blup_py(zi,xi,k)=v2
        end do !k
      end do !zi
    end do !xi

    !maternal environmental effect (mee)
    meen=0
    do zi=1,pedn   !same as rn
      if (dsgn(zi,2) .ne. 0) then
        do j=1,meen    !maternal environment effect #
          if (dsgn(zi,2) == mee(j)) goto 77
        end do
        meen=meen+1
        mee(meen)=dsgn(zi,2)
      end if
77    continue
    end do
    !print*,meen
    !BLUP for mee
    do xi=1,tn
      do zi=1,meen
        v1=0
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          do j=1,rnm(yi)
            if (dsgn2(yi,j,2)==mee(zi)) then
              v1=v1+vcva(4,xi,yi)*py(trny-rnm(yi)+j,1)
            end if
          end do
        end do
        blup_mee(zi,xi)=v1
      end do
    end do

  !open (unit=45,file=trim(fl11)//".py",status='unknown')
  !do xi=1,tn
  !  do i=1,pedn
  !    write(45,'(i3,100f24.16)')xi,blup_py(i,xi,1:2)
  !  end do
  !end do
  !close(45)

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
      write(43,'(i3,100f24.16)')xi,blup_ebv(i,xi,1:2)
    end do
  end do
  close(43)
  if (mn==4) then
    write(43,*)''
    write(43,'(a5,a26)')'trait','MEE'
    do xi=1,tn
      do i=1,meen
        write(43,'(i3,100f24.16)')xi,blup_mee(i,xi)
      end do
    end do
  end if
  close(43)

end if


if (fl11r.ne.'null') then
  print*,'for BLUP solution and reliability after convergence ***********'
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
        !do k=1,mn
        do k=1,2
          v1=0;v2=0
          trny=0
          do yi=1,tn
            trny=trny+rnm(yi)
            do j=1,rnm(yi)
              if (dsgn2(yi,j,k).ne.0) then
                v1=v1+mbin(1,dsgn2(yi,j,k),zi)*vcva(k,xi,yi)*py(trny-rnm(yi)+j,1)
              end if
              if (dsgn2(yi,j,(3-k)).ne.0) then
                v1=v1+mbin(1,dsgn2(yi,j,(3-k)),zi)*vcva(3,xi,yi)*py(trny-rnm(yi)+j,1)
              end if
              if (dsgn2(yi,j,k)==zi) then
                v2=v2+vcva(k,xi,yi)*py(trny-rnm(yi)+j,1)  !not sure for this >>> check
              end if
            end do
          end do
          blup_ebv(zi,xi,k)=v1
          blup_py(zi,xi,k)=v2
        end do !k
      end do !zi
    end do !xi

    !to get reliability for BLUP, i.e. GPA
    do vi=1,tn
      !$OMP PARALLEL DO PRIVATE(zi,k,v1,v2)
      do zi=1,pedn
        !do k=1,mn
        do k=1,2
          v2=0;v22=0;v222=0
          trnx=0
          do xi=1,tn
            trnx=trnx+rnm(xi)
            do i=1,rnm(xi)
              v1=0;v11=0
              trny=0
              do yi=1,tn
                trny=trny+rnm(yi)
                do j=1,rnm(yi)
                  if (dsgn2(yi,j,k).ne.0) then
                    v1=v1+mbin(1,dsgn2(yi,j,k),zi)*vcva(k,vi,yi)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
                  end if
                  if (dsgn2(yi,j,(3-k)).ne.0) then
                    v11=v11+mbin(1,dsgn2(yi,j,(3-k)),zi)*vcva(3,vi,yi)*pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)
                  end if
                  !print*,v1,v11,dsgn(j,k),i,j,k,vi,yi,zi
                end do
              end do
              if (dsgn2(xi,i,k).ne.0) then
                v2=v2+v1*mbin(1,dsgn2(xi,i,k),zi)*vcva(k,vi,xi)
              end if
              if (dsgn2(xi,i,(3-k)).ne.0) then
                v22=v22+v11*mbin(1,dsgn2(xi,i,(3-k)),zi)*vcva(3,vi,xi)
                v222=v222+2*(v1*mbin(1,dsgn2(xi,i,(3-k)),zi)*vcva(3,vi,xi))
              end if
              !print*,v2,v22,v222,dsgn(i,k),zi,rnm(xi),tn,rnm(yi)
              !pause
            end do
          end do
          !blup_r(zi,vi,k)=v2
          blup_r(zi,vi,k)=(v2+v22+v222)/vcva(k,vi,vi)
          !print*,v2,v22,v222
          !pause
        end do
      end do
      !$OMP END PARALLEL DO 
    end do

    !maternal environmental effect (mee)
    meen=0
    do zi=1,pedn   !same as rn
      if (dsgn(zi,2) .ne. 0) then
        do j=1,meen    !maternal environment effect #
          if (dsgn(zi,2) == mee(j)) goto 78
        end do
        meen=meen+1
        mee(meen)=dsgn(zi,2)
      end if
78    continue
    end do
    !print*,meen
    !BLUP for mee
    do xi=1,tn
      do zi=1,meen
        v1=0
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          do j=1,rnm(yi)
            if (dsgn2(yi,j,2)==mee(zi)) then
              v1=v1+vcva(4,xi,yi)*py(trny-rnm(yi)+j,1)
            end if
          end do
        end do
        blup_mee(zi,xi)=v1
      end do
    end do
  do vi=1,tn
    !$OMP PARALLEL DO PRIVATE(zi,k,v1,v2)
    do zi=1,meen
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
              if (dsgn2(yi,j,2)==mee(zi)) then
                v1=v1+vcva(4,vi,yi)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
              end if
            end do
          end do
          if (dsgn2(xi,i,2)==mee(zi)) then
            !v2=v2+v1
            v2=v2+v1*vcva(4,vi,xi)
          end if
        end do
      end do
      !blup_mee_r(zi,vi)=v2
      blup_mee_r(zi,vi)=v2/vcva(4,vi,vi)
    end do
    !$OMP END PARALLEL DO 
  end do

  !open (unit=45,file=trim(fl11r)//".py",status='unknown')
  !do xi=1,tn
  !  do i=1,pedn
  !    write(45,'(i3,100f24.16)')xi,blup_py(i,xi,1:2)
  !  end do
  !end do
  !close(45)

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
      write(43,'(i3,100f24.16)')xi,blup_ebv(i,xi,1:2)
    end do
  end do
  if (mn==4) then
    write(43,*)''
    write(43,'(a5,a26)')'trait','MEE'
    do xi=1,tn
      do i=1,meen
        write(43,'(i3,100f24.16)')xi,blup_mee(i,xi)
      end do
    end do
  end if
  close(43)

  open (unit=46,file=trim(fl11r)//".r2",status='unknown')
  do xi=1,tn
    do i=1,pedn
      !write(46,'(i3,100f24.16)')xi,blup_r(i,xi,:)**2
      write(46,'(i3,100f24.16)')xi,blup_r(i,xi,1:2)
    end do
  end do
  if (mn==4) then
    write(46,*)''
    do xi=1,tn
      do i=1,meen
        write(46,'(i3,100f24.16)')xi,blup_mee_r(i,xi)
      end do
    end do
  end if
  close(46)

end if


end subroutine

!end module
