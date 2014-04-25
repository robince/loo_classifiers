!
! Linear (nearest mean) decoding with leave one out cross validation
! Cesare v2 algorithm
!
#include "fintrf.h"

#define DATATYPE c_float
#define MATDATATYPE mxSINGLE_CLASS

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use MatlabAPImex
    use MatlabAPImx
    use iso_c_binding
    implicit none
! ARG
    mwPointer plhs(*), prhs(*)
    integer(4) nlhs, nrhs
! LOC
    real(kind=DATATYPE), pointer :: X(:,:), conmtx(:,:)
    integer(kind=c_int16_t), pointer :: Y(:), prdY(:)
    real(kind=DATATYPE), allocatable :: sum_data(:,:), curtrl(:), stmdif(:,:), stmdst(:)
    real(kind=DATATYPE), allocatable :: Ntrlcls(:), prctrl(:), trlratio(:)
    real(kind=c_double), pointer :: Ncls_r
    mwSize :: thscls
    mwSize, allocatable :: minidx(:)
    mwSize, parameter :: One=1
    mwSize :: Nftr, Ncls, Ntrl, fi, ti, ssi

    ! Process Inputs
    if( nrhs /= 3 ) then
        call mexErrMsgTxt("This function takes 3 inputs")
    endif
    if( nlhs /= 2) then
        call mexErrMsgTxt("This function returns 2 outputs")
    endif

    call fpGetPr( X, prhs(1) )
    call fpGetPr( Y, prhs(2) )
    call fpGetPr( Ncls_r, prhs(3) )

    if( (.not.associated(X)) .or. (.not.associated(Y)) .or. (.not.associated(Ncls_r)) ) then
        call mexErrMsgTxt("Problem with inputs: check types and dimensions")
    endif
    Ncls = Ncls_r

    Nftr = size(X,1)
    Ntrl = size(X,2)

    plhs(1) = mxCreateNumericMatrix(Ncls, Ncls, MATDATATYPE, mxREAL)
    call fpGetPr( conmtx, plhs(1) )

    plhs(2) = mxCreateNumericMatrix(Ntrl, One, mxINT16_CLASS, mxREAL)
    call fpGetPr( prdY, plhs(2) )
    
    conmtx = 0.0

    allocate(sum_data(Nftr, Ncls))
    allocate(Ntrlcls(Ncls))
    Ntrlcls = 0.0;
    sum_data = 0.0;

    ! sum data
    do ti=1,Ntrl
      Ntrlcls(Y(ti)) = Ntrlcls(Y(ti)) + 1
      do fi=1,Nftr
        sum_data(fi, Y(ti)) = sum_data(fi, Y(ti)) + X(fi,ti)
      end do
    end do

    allocate(curtrl(Nftr))
    allocate(stmdif(Nftr,Ncls))
    allocate(stmdst(Ncls))
    allocate(minidx(1))
    allocate(prctrl(Nftr))
    allocate(trlratio(Nftr))

    prctrl = 100.0 / Ntrlcls
    trlratio = Ntrlcls / (Ntrlcls - 1.0)

    do ti=1,Ntrl
      thscls = Y(ti);
      curtrl = X(:,ti) * Ntrlcls(thscls)
      do ssi=1,Ncls
          do fi=1,Nftr
              stmdif(fi,ssi) = sum_data(fi,ssi) - curtrl(fi)
          end do
      end do
      stmdif(:,thscls) = stmdif(:,thscls) * trlratio(thscls)
      do ssi=1,Ncls
          stmdst(ssi) = sum(stmdif(:,ssi) * stmdif(:,ssi))
      end do
      minidx = minloc(stmdst)
      prdY(ti) = minidx(1)
      conmtx(thscls, minidx(1)) = conmtx(thscls, minidx(1)) + prctrl(thscls)
    end do

    deallocate(sum_data)
    deallocate(Ntrlcls)
    deallocate(curtrl)
    deallocate(stmdif)
    deallocate(stmdst)
    deallocate(minidx)
    deallocate(prctrl)
    deallocate(trlratio)

end subroutine mexFunction

