!
! Diagonal Linear
!
#include "fintrf.h"


subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use MatlabAPImex
    use MatlabAPImx
    implicit none
! ARG
    mwPointer plhs(*), prhs(*)
    integer(4) nlhs, nrhs
! LOC
    real(8), pointer :: dat(:,:,:), conmtx(:,:), prdstm(:,:)
    real(8), allocatable :: ftrsum(:,:), curtrl(:), stmdif(:,:), stmdst(:)
    real(8), allocatable :: ftravg(:,:), curftravg(:,:), u(:), v(:)
    real(8), allocatable :: cendat(:,:,:), xc(:,:), ftrvar(:), curftrvar(:)
    real(8) :: prctrl, trlfac
    mwSize, allocatable :: minidx(:)
    mwSize :: Nftr, Ncls, Ntrl, fi, ci, ti, cci, Ntottrl, Ntrl1, dof

    ! Process Inputs
    if( nrhs /= 1 ) then
        call mexErrMsgTxt("This function takes 1 inputs")
    endif
    if( nlhs /= 2) then
        call mexErrMsgTxt("This function returns 2 output")
    endif

    call fpGetPr( dat, prhs(1) )

    if( (.not.associated(dat)) ) then
        call mexErrMsgTxt("Problem with inputs: check types and dimensions")
    endif

    Nftr = size(dat,1)
    Ncls = size(dat,2)
    Ntrl = size(dat,3)
    Ntottrl = Ntrl * Ncls
    dof = Ntottrl - 1 - Ncls
    Ntrl1 = Ntrl - 1

    plhs(1) = mxCreateNumericMatrix( Ncls, Ncls, mxDOUBLE_CLASS, mxREAL)
    call fpGetPr( conmtx, plhs(1) )

    plhs(2) = mxCreateNumericMatrix( Ntrl, Ncls, mxDOUBLE_CLASS, mxREAL)
    call fpGetPr( prdstm, plhs(2) )

    conmtx = 0.0d0

    allocate(ftrsum(Nftr, Ncls))
    allocate(ftravg(Nftr, Ncls))
    allocate(ftrvar(Nftr))
    allocate(curftrvar(Nftr))
    allocate(curftravg(Nftr, Ncls))
    allocate(curtrl(Nftr))
    allocate(u(Nftr))
    allocate(v(Nftr))
    allocate(stmdif(Nftr,Ncls))
    allocate(stmdst(Ncls))
    allocate(minidx(1))
    allocate(cendat(Nftr,Ncls,Ntrl))
    allocate(xc(Nftr, Ncls*Ntrl))

    ftrsum = sum(dat,3)
    ftravg = ftrsum / Ntrl
    do ti=1,Ntrl
        cendat(:,:,ti) = dat(:,:,ti) - ftravg
    end do
    xc = reshape(cendat, (/ Nftr, Ncls*Ntrl /) )
    ftrvar = sum(xc*xc, 2)
    prctrl = 100.0d0 / Ntrl

    do ci=1,Ncls
        curftravg = ftravg
        do ti=1,Ntrl
            curtrl = dat(:,ci,ti)
            curftravg(:,ci) = (ftrsum(:,ci) - curtrl) / Ntrl1
            u = curtrl - ftravg(:,ci)
            v = curtrl - curftravg(:,ci)
            curftrvar = (ftrvar - (u*v)) / dof

            do cci=1,Ncls
                u = curtrl - curftravg(:,cci)
                stmdst(cci) = sum((u*u) / curftrvar)
            end do

            minidx = minloc(stmdst)
            prdstm(ti,ci) = minidx(1)

            conmtx(ci, minidx(1)) = conmtx(ci, minidx(1)) + prctrl
        end do
    end do

    deallocate(ftrsum)
    deallocate(ftravg)
    deallocate(ftrvar)
    deallocate(curftrvar)
    deallocate(curftravg)
    deallocate(curtrl)
    deallocate(u)
    deallocate(v)
    deallocate(stmdif)
    deallocate(stmdst)
    deallocate(minidx)
    deallocate(cendat)

end subroutine mexFunction
