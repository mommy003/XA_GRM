
subroutine submat (mbin,yidx,yv,rn,rnm,trn,pedn,tn,fnm,tfn,mn,subgrm)        

!***********************************************************************
!subgrm: sub matrix extraction  
!S. Hong Lee (2021)

!mbin : matrices bin
!obs  : observed ID
!yv   : phenotypes
!rn   : no. phenotypes
!pedn : no. pedigree (ID) should be = rn
!fn   : no. fixed effects
!mn   : no. matrices (random effects)
!***********************************************************************

implicit none

INTEGER::n,nped,ix,iz,tn,nb,nbc,nbc2,xi,yi,ri,zi,vi,ui,wi,wj,nit,i,j
integer::rn,pedn,mn,trn,trnx,trny,tfn,trnv,trnu,itit
integer::yidx(tn,rn),rnm(tn),fnm(tn)
real::mbin(mn,pedn,pedn),yv(rn,tn),sv_fact
double precision::tyv(trn,1)

character(len=64)::subgrm

!double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),conv
!double precision::vcva(mn,tn,tn),vcve(tn,tn),vvc(tn)
!double precision::fix_vcva(mn,tn,tn),fix_vcve(tn,tn)
!double precision::vcva_rsv(mn,tn,tn),vcve_rsv(tn,tn)
!double precision::blup_ebv(pedn,tn,mn),blup_py(pedn,tn,mn),beta(tfn,1) !,cor
!double precision::blup_r(pedn,tn,mn)
!double precision,allocatable::tmp_v(:,:),gpg(:,:),grm(:,:)

!integer::i,j,m,k,l,io
!double precision::sum_h2,sum_v(tn)


!double precision::LC,LA,LG,ypy,x1,x2,x3,x10,x11,v1,v2,v10,v11
!double precision::y1,y2,y3,y10,y11,z1,z2,z3,z10,z11
!double precision::LKHP,MLKH,mva(mn),mva2(mn),mve,mve2,mcov(mn)
!double precision::h2(mn),h2n(mn),tr1,tr2(mn),tr3

double precision::xm(trn,tfn),xmm(tn,pedn,10000)

!double precision::xvmx(tfn,tfn),xvm(tfn,trn) 
!double precision::xvm2(tfn,trn),xvm3(trn,tfn) 
!double precision::pm(trn,trn),py(trn,1),xb(trn,1)

!double precision::v_tmp(trn),v2_tmp(trn)

!double precision::aim(tn+mn*(tn+(tn**2-tn)/2),tn+mn*(tn+(tn**2-tn)/2))
!double precision::sdmI(tn+mn*(tn+(tn**2-tn)/2),tn+mn*(tn+(tn**2-tn)/2))
!double precision::up(tn+mn*(tn+(tn**2-tn)/2),1),dldv(tn+mn*(tn+(tn**2-tn)/2),1)

!for BLUP ebv
!double precision::vmat(mn,pedn,rn)

!!time measure **********************************************************
!INTEGER::now(8)
!CHARACTER*8 date
!CHARACTER*20 time,zone
!double precision:: timef, t1_real, t2_real, elapsed_real_secs
!double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
!!***********************************************************************

!character(len=7)::cdum
!character(len=64)::fl4,fl24,fl5,fl11,fl11r,fl11rr,fl11rrr,wt_res,cdum_t(tn)
!logical::file_exist

!!CDF
!real ( kind = 8 ) bound,mean,p,q,sd,x
!integer ( kind = 4 ) status,which

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

  open (UNIT=41,FILE=subgrm,STATUS='unknown')
  do i=1,rnm(1)
    do j=1,i
      !print*,yidx(1,i),yidx(1,j),mbin(1,yidx(1,i),yidx(1,j))
      write(41,'(i8,i8,f14.8)')i,j,real(mbin(1,yidx(1,i),yidx(1,j)))
    end do
    !pause
  end do
  close(41)
  print*,'Sub matrix (for non-missing phenotyped IDs) stored in ',trim(subgrm)
  print*,'The order follows in -d file. Also see ',trim(subgrm)//".dat"

  open (UNIT=42,FILE=trim(subgrm)//".dat",STATUS='unknown')
  do i=1,rnm(1)
    write(42,'(i8,i8,f24.16)')i,i,tyv(i,1)
  end do
  close(42)



end subroutine

