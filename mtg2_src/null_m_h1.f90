
subroutine null_m_h (yidx,yv,rn,rnm,trn,pedn,tn,fnm,tfn,mn,fl4,fl5,fl11,xmm,nit,conv)        

!***********************************************************************
!aireml_b: bivariate analysis especially for case control 
!S. Hong Lee

!mbin : matrices bin
!obs  : observed ID
!yv   : phenotypes
!rn   : no. phenotypes
!pedn : no. pedigree (ID) should be = rn
!fn   : no. fixed effects should be = 1 (mean)
!mn   : no. matrices (random effects)
!***********************************************************************

implicit none

INTEGER::n,nped,ix,iz,tn,nb,nbc,nbc2,xi,yi,ri,zi,vi,ui,wi,wj,nit
integer::rn,pedn,mn,trn,trnx,trny,tfn,trnv,trnu
integer::yidx(tn,rn),rnm(tn),fnm(tn)
integer::rn1,rn2,yidx1(rn),yidx2(rn)
!double precision::mbin(mn,pedn,pedn),yv(rn,2),tyv(rn*2,1)
real::mbin(mn,pedn,pedn),yv(rn,tn),tyv(trn,1)

!double precision::qq(mn,rn,rn),qq2(mn,rn,rn)
!real::qq(mn,rn,rn),qq2(mn,rn,rn)

double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),conv
double precision::vcva(mn,tn,tn),vcve(tn,tn),vvc(tn)
double precision::blup_ebv(pedn,1),beta(tfn,1) !,cor
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
character(len=64)::fl4,fl5,fl11
logical::file_exist


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

print*,'*** number of records used ***'
do i=1,tn
  print*,'trait',i,':',rnm(i)
end do
print*,''



  ! y = Xb + Zu + e ---- fitting effects without QTL

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
    vcve(xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))
    !vcva(:,xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))/(2*mn)
  end do

  ! Start iteration
  LKH=-1000000000
  MLKH=-1000000000


  ! V matrix    V= ZAZ'Va+Z2GZ2'Vpe+IVe*************************
  !do zi=1,20
  do zi=1,nit
    LKHP=LKH
  !LKH and AIREML (without sparce technique)
  ! V matrix    V = ZAZ'Va + IVe *************************

  call cpu_time (t2_cpu)

  pm=0;trnx=0!;LA=0
  !do xi=1,tn
  !  trnx=trnx+rnm(xi)
  !  do i=1,rnm(xi)
  !    pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)=1/vcve(xi,xi)
  !    LA=LA+log(vcve(xi,xi))
  !  end do
  !end do
  do xi=1,tn
    trnx=trnx+rnm(xi)
    trny=0
    do yi=1,xi
      trny=trny+rnm(yi)
      !$OMP PARALLEL DO PRIVATE(i,j,k) 
      do i=1,rnm(xi)
        do j=1,rnm(yi)
          if (yidx(xi,i)==yidx(yi,j)) then
            !diagonal for the same ID, off-diagonal (cov) for the same ID
            pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)+vcve(xi,yi)
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
  !print*,LA

    ! P = (VI)-(VI X (X'VI X)I' X' VI)******************************
    !xvm=MATMUL(transpose(xm),pm)
    do i=1,tfn
      do j=1,trn
        v1=0
        do k=1,trn
          v1=v1+xm(k,i)*pm(k,j)    !xm*pm
        end do
        xvm(i,j)=v1
      end do
    end do

    !xvmx=MATMUL(xvm,xm)
    do i=1,tfn
      do j=1,tfn
        v1=0
        do k=1,trn
          v1=v1+xvm(i,k)*xm(k,j)    !xm*pm
        end do
        xvmx(i,j)=v1
      end do
    end do

    !call cholesky_inv (xvmx,tfn,xvmxI,LG)
    call cholesky_inv (xvmx,tfn,LG)
    !BLUE = (X'VI X)I XVI y **********************
    !beta=matmul(matmul(matmul(xvmx,transpose(xm)),pm),tyv)
    do i=1,tfn
      do j=1,trn
        v1=0
        do k=1,tfn
          v1=v1+xvmx(i,k)*xm(j,k)    !xvmx*xm
        end do
        xvm2(i,j)=v1
      end do
    end do

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
   !*********************************************
   !print*,'beta:',beta

    !vmx4=MATMUL(MATMUL(MATMUL(vmI,xm),xvmxI),MATMUL(transpose(xm),vmI))
    !pm=MATMUL(MATMUL(MATMUL(vmI,xm),xvmxI),MATMUL(transpose(xm),vmI))
    !do i=1,trn
    !  do j=1,trn
    !    !pm(i,j)=vmI(i,j)-vmx4(i,j)
    !    !pm(i,j)=vmI(i,j)-pm(i,j)
    !    pm(i,j)=vm(i,j)-pm(i,j)
    !  end do
    !end do
    !pm=pm-MATMUL(MATMUL(MATMUL(pm,xm),xvmx),MATMUL(transpose(xm),pm))
    do i=1,trn
      do j=1,trn
        v1=0
        do k=1,tfn
          v1=v1+xvm3(i,k)*xvm(k,j)
        end do
        pm(i,j)=pm(i,j)-v1
      end do
    end do

    !ypy estimation
    !ypyf=MATMUL(MATMUL(transpose(tyv),pm),tyv)
    v2=0
    do i=1,trn
      v1=0
      do j=1,trn
        v1=v1+tyv(j,1)*pm(j,i)
      end do
      v2=v2+v1*tyv(i,1)
    end do

    !LKH=-0.5*(LA+LG+ypyf(1,1))
    LKH=-0.5*(LA+LG+v2)
    !print*,LA,LG,v2  !ypyf

    call cpu_time (t2_cpu)
    !print*,'LKH:',t2_cpu-t1_cpu

    PRINT '(a7,100f12.4)','LKH',LKH
    do xi=1,tn
      PRINT '(a7,100f12.4)','Ve',vcve(xi,xi)
    end do
    print*,''
    do wi=1,mn
      do xi=1,tn
        PRINT '(a7,100f12.4)','Va',vcva(wi,xi,xi)
      end do
      do xi=2,tn
        do yi=1,xi-1
          PRINT '(a7,100f12.4)','cov',vcva(wi,xi,yi)
        end do
      end do
      print*,''
    end do
    !pause


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
          v1=v1+py(trnx-rnm(xi)+j,1)*pm(trnx-rnm(xi)+j,trnx-rnm(xi)+i)
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
        tr1=tr1+pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)
      end do
      tr1=-0.5*tr1
      v2=0
      do i=1,rnm(xi)
        v2=v2+py(trnx-rnm(xi)+i,1)*py(trnx-rnm(xi)+i,1)
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
            v1=v1+py(trny-rnm(yi)+j,1)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
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
            v2=v2+v2_tmp(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)                                !*py
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
              v2=v2+v2_tmp(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)                                !*py
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
    print*,'derivatives done - time:',real(t2_cpu-t1_cpu)
    print*,''
 
    !do i=1,(2+mn*3)
    do i=1,tn+nb*mn
      aim(i,i)=aim(i,i)/2
      do j=1,i-1
        aim(i,j)=aim(i,j)/2
        aim(j,i)=aim(i,j)
      end do
    end do

    !do i=1,tn+nb*mn
    !  print*,real(aim(i,1:i))*2
    !end do
    !pause

    !call cholesky_inv (aim,(2+mn*3),x1)
    call cholesky_inv (aim,(tn+nb*mn),x1)
    print *,''

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
    !do i=1,tn+nb*mn
    !  print*,real(up(i,1)),real(dldv(i,1))
    !end do
    !print*,''
    !pause

    !ve=ve+up(1,1)
    !ve2=ve2+up(2,1)
    !do xi=1,mn
    !  va(xi)=va(xi)+up(2+xi*3-2,1)
    !  va2(xi)=va2(xi)+up(2+xi*3-1,1)
    !  cov(xi)=cov(xi)+up(2+xi*3,1)
    !end do
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

    !print*,ve,ve2,va,va2,cov
    !print*,vcve,vcva

  END do !zi



  !PRINT*,'LKH not converged na na na'
  !WRITE(41,*)'LKH not converged na na na'
  va=mva;va2=mva2;ve=mve;ve2=mve2;LKH=MLKH;cov=mcov


1000 continue 


  open (UNIT=41,FILE=fl5,STATUS='unknown')

  do xi=1,tn
    write(41,'(a7,100f12.4)')'Ve',vcve(xi,xi),sqrt(sdmI(xi,xi))
  end do
  write(41,*)''

  do xi=1,tn
    sum_v(xi)=vcve(xi,xi)
  end do
  do wi=1,mn
    do xi=1,tn
      write(41,'(a7,100f12.4)')'Va',vcva(wi,xi,xi),sqrt(sdmI(tn+wi*nb-nb+xi,tn+wi*nb-nb+xi))
      sum_v(xi)=sum_v(xi)+vcva(wi,xi,xi)
    end do
    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        write(41,'(a7,100f12.4)')'cov',vcva(wi,xi,yi),sqrt(sdmI(tn+wi*nb-nb+tn+nbc,tn+wi*nb-nb+tn+nbc))
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
      !x1=va(i)/sum_v
      x1=vcva(wi,xi,xi)/sum_v(xi)
      !x11=sdmI(2+i*3-2,1) !for the row ************
      x11=sdmI(tn+nb*wi-nb+xi,xi) !for the row ************

      do wj=1,mn
        !x11=x11+sdmI(2+i*3-2,2+j*3-2)
        x11=x11+sdmI(tn+nb*wi-nb+xi,tn+nb*wj-nb+xi)
      end do
      !print*,x11

      !SE **********************************
      !x2=x1**2*(sdmI(2+i*3-2,2+i*3-2)/va(i)**2+x10/sum_v**2-2*x11/(va(i)*sum_v))
      x2=x1**2*(sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+xi)/vcva(wi,xi,xi)**2+vvc(xi)/sum_v(xi)**2-2*x11/(vcva(wi,xi,xi)*sum_v(xi)))
      write(41,'(a7,100f12.4)')'h2',x1,sqrt(x2)

    end do
    write(41,*)''

    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        !z1=cov(i)/sqrt(va(i)*va2(i))
        z1=vcva(wi,xi,yi)/sqrt(vcva(wi,xi,xi)*vcva(wi,yi,yi))
        !var(cor) *************************************
        !z2=sdmI(2+i*3-2,2+i*3-2)/(4*va(i)**2)
        z2=sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+xi)/(4*vcva(wi,xi,xi)**2)
        !z2=z2+sdmI(2+i*3-1,2+i*3-1)/(4*va2(i)**2)
        z2=z2+sdmI(tn+nb*wi-nb+yi,tn+nb*wi-nb+yi)/(4*vcva(wi,yi,yi)**2)
        !z2=z2+sdmI(2+i*3,2+i*3)/(cov(i)**2)
        z2=z2+sdmI(tn+nb*wi-nb+tn+nbc,tn+nb*wi-nb+tn+nbc)/(vcva(wi,xi,yi)**2)
        !z2=z2+2*sdmI(2+i*3-2,2+i*3-1)/(4*va(i)*va2(i))
        z2=z2+2*sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+yi)/(4*vcva(wi,xi,xi)*vcva(wi,yi,yi))
        !z2=z2-2*sdmI(2+i*3-2,2+i*3)/(2*va(i)*cov(i))
        z2=z2-2*sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+tn+nbc)/(2*vcva(wi,xi,xi)*vcva(wi,xi,yi))
        !z2=z2-2*sdmI(2+i*3,2+i*3-1)/(2*cov(i)*va2(i))
        z2=z2-2*sdmI(tn+nb*wi-nb+tn+nbc,tn+nb*wi-nb+yi)/(2*vcva(wi,xi,yi)*vcva(wi,yi,yi))
        z2=z2*(z1**2)
        write(41,'(a7,100f12.4)')'cor',z1,sqrt(z2)
      end do
    end do
    write(41,*)''
  end do

  write(41,*)''
  write(41,'(a7,100f12.4)')'LKH',LKH

  write(41,*)''
  do i=1,(tn+nb*mn)
    write(41,*)real(sdmI(i,1:i))
  end do

  close(41)


  xb=tyv-matmul(xm,beta)
  !open (unit=43,file='ascm.res',status='unknown')
  open (unit=43,file=trim(fl11)//".res",status='unknown')
  k=0
  do xi=1,tn
    do i=1,pedn
      if (yv(i,xi).ne.-99999) then
        k=k+1
        write(43,*)xi,xb(k,1)
      else
        write(43,*)xi,'NA'
      end if
    end do
  end do
  close(43)


if (fl11.ne.'null') then

  print*,'for BLUP solution after convergence *******************'

  open (unit=47,file=trim(fl11)//".fsl",status='unknown')
  do xi=1,tfn
    write(47,'(i3,100f24.16)')xi,beta(xi,1),sqrt(xvmx(xi,xi))
  end do
  close(47)

end if



end subroutine

