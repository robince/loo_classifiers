!
! Linear (nearest mean) decoding with leave one out cross validation
! Cesare v2 algorithm
!
#include "fintrf.h"


subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use MatlabAPImex
    use MatlabAPImx
    !use lapack95, only: potrf, potri
    !use blas95, only: gemv, gemm, gem2vu, dot
    use lapack95
    use blas95
    implicit none
! ARG
    mwPointer plhs(*), prhs(*), lhs(1), rhs(1)
    integer(4) nlhs, nrhs
! LOC
    real(8), pointer :: dat(:,:,:), conmtx(:,:), prdstm(:,:), invftrcov(:,:)
    real(8), allocatable :: ftravg(:,:), ftrsum(:,:),curtrl(:), stmdif(:,:), stmdst(:)
    real(8), allocatable :: cendat(:,:,:), ftrcov(:,:), xc(:,:)
    real(8), allocatable :: facinvftrcov(:,:), curftravg(:,:), curinvftrcov(:,:)
    real(8), allocatable :: u(:), v(:), Au(:), vtA(:), num(:,:)
    real(8) :: prctrl, invfac, den
    mwSize, allocatable :: minidx(:)
    mwSize :: Nftr, Ncls, Ntrl, Ntrl1, Ntottrl, Ntottrl1, fi, ci, ti, cci, k

    ! Process Inputs
    if( nrhs /= 1 ) then
        call mexErrMsgTxt("This function takes 1 input")
    endif
    if( nlhs /= 2) then
        call mexErrMsgTxt("This function returns 2 output")
    endif

    dat => fpGetPr3( prhs(1) )

    if( (.not.associated(dat)) ) then
        call mexErrMsgTxt("Problem with inputs: check types and dimensions")
    endif

    Nftr = size(dat,1)
    Ncls = size(dat,2)
    Ntrl = size(dat,3)
    Ntrl1 = Ntrl - 1
    Ntottrl = Ntrl * Ncls
    Ntottrl1 = Ntottrl - 1

    conmtx => fpAllocate2( Ncls, Ncls )
    prdstm => fpAllocate2( Ntrl, Ncls )
    conmtx = 0.0d0

    plhs(1) = mxArrayHeader(conmtx)
    plhs(2) = mxArrayHeader(prdstm)
    
    ! mean responses
    allocate(ftravg(Nftr, Ncls))
    allocate(curftravg(Nftr, Ncls))
    allocate(ftrsum(Nftr, Ncls))
    ftrsum = sum(dat,3)
    ftravg = ftrsum / Ntrl;
    ! center data
    allocate(cendat(Nftr,Ncls,Ntrl))
    do ti=1,Ntrl
        cendat(:,:,ti) = dat(:,:,ti) - ftravg
    end do

    allocate(ftrcov(Nftr, Nftr))
    allocate(facinvftrcov(Nftr, Nftr))
    allocate(xc(Nftr, Ncls*Ntrl))
    xc = reshape(cendat, (/ Nftr, Ncls*Ntrl /) )
    ftrcov = matmul(xc, transpose(xc)) / Ntrl1

    ! would like to use mxArrayHeader(Destroy) here
    ! to avoid copy, but get crashes when I do that
    ! and don't have time to investigate
    rhs(1) = mxArray(ftrcov)
    k = mexCallMATLAB(1, lhs, 1, rhs, "inv")
    call mxDestroyArray(rhs(1))
    invftrcov => fpGetPr(lhs(1))

    !allocate(curtrl(Nftr))
    allocate(stmdif(Nftr,Ncls))
    allocate(stmdst(Ncls))
    allocate(minidx(1))
    allocate(u(Nftr))
    allocate(v(Nftr))
    allocate(Au(Nftr))
    allocate(vtA(Nftr))
    allocate(num(Nftr,Nftr))
    allocate(curinvftrcov(Nftr,Nftr))

    prctrl = 100.0d0 / Ntrl
    invfac = (Ntottrl - 2.0d0) / (Ntottrl - 1.0d0)

    facinvftrcov = invfac * invftrcov
    do ci=1,Ncls
        curftravg = ftravg
        do ti=1,Ntrl
            !curtrl = dat(:,ci,ti)
            ! update mean
            curftravg(:,ci) = (ftrsum(:,ci)-dat(:,ci,ti)) / Ntrl1
            u = ftravg(:,ci) - dat(:,ci,ti)
            v = dat(:,ci,ti) - curftravg(:,ci)

            ! A*u
            !call gemv(invftrcov, u, Au, 1.0d0, 0.0d0, 'N') 
            ! v'*A
            !call gemv(invftrcov, v, vtA, 1.0d0, 0.0d0, 'T') 
            call gem2v(invftrcov, u, v, Au, vtA, 1.0d0, 0.0d0)
            ! num = (A*u)*(v'*A)
            num = spread(Au,dim=2,ncopies=Nftr)*spread(vtA,dim=1,ncopies=Nftr)

            ! reuse Au variable for Av
            call gemv(invftrcov, v, Au, 1.0d0, 0.0d0, 'N')
            den = (Ntottrl1 + dot(u, Au)) / invfac
            curinvftrcov = facinvftrcov - (num / den)

            do cci=1,Ncls
                u = dat(:,ci,ti) - curftravg(:,cci)
                call gemv(curinvftrcov, u, Au, 1.0d0, 0.0d0, 'N')
                stmdst(cci) = dot(u, Au)
            end do

            minidx = minloc(stmdst)
            prdstm(ti,ci) = minidx(1)

            conmtx(ci, minidx(1)) = conmtx(ci, minidx(1)) + prctrl
            
        end do
    end do

    deallocate(u)
    deallocate(v)
    deallocate(Au)
    deallocate(vtA)
    deallocate(num)
    deallocate(xc)
    deallocate(ftrcov)
    deallocate(curinvftrcov)
    deallocate(facinvftrcov)
    deallocate(curftravg)
    deallocate(ftravg)
    deallocate(cendat)
    deallocate(ftrsum)
    !deallocate(curtrl)
    deallocate(stmdif)
    deallocate(stmdst)
    deallocate(minidx)
    nullify(invftrcov)
    call mxDestroyArray(lhs(1))


end subroutine mexFunction

