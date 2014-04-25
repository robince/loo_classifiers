!
! Linear (nearest mean) decoding with leave one out cross validation
! Cesare v2 algorithm
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
    real(8), pointer :: resp(:,:,:), conmtx(:,:), prdstm(:,:)
    real(8), allocatable :: sumresp(:,:), currsp(:), stmdif(:,:), stmdst(:)
    real(8) :: prctrl, trlfac
    mwSize, allocatable :: minidx(:)
    mwSize :: Nftr, Nstm, Ntrl, fi, si, ti, ssi

    ! Process Inputs
    if( nrhs /= 1 ) then
        call mexErrMsgTxt("This function takes 2 inputs")
    endif
    if( nlhs /= 2) then
        call mexErrMsgTxt("This function returns 2 output")
    endif

    resp => fpGetPr3( prhs(1) )

    if( (.not.associated(resp)) ) then
        call mexErrMsgTxt("Problem with inputs: check types and dimensions")
    endif

    Nftr = size(resp,1)
    Nstm = size(resp,2)
    Ntrl = size(resp,3)

    conmtx => fpAllocate2( Nstm, Nstm )
    prdstm => fpAllocate2( Ntrl, Nstm )

    plhs(1) = mxArrayHeader(conmtx)
    plhs(2) = mxArrayHeader(prdstm)
    
    conmtx = 0.0d0

    allocate(sumresp(Nftr, Nstm))
    allocate(currsp(Nftr))
    allocate(stmdif(Nftr,Nstm))
    allocate(stmdst(Nstm))
    allocate(minidx(1))

    sumresp = sum(resp,3)
    prctrl = 100.0d0 / Ntrl
    trlfac = real(Ntrl,8) / (Ntrl - 1.0d0)

    do si=1,Nstm
        do ti=1,Ntrl
            currsp = Ntrl * resp(:,si,ti)
            do ssi=1,Nstm
                do fi=1,Nftr
                    stmdif(fi,ssi) = sumresp(fi,ssi) - currsp(fi)
                end do
            end do
            stmdif(:,si) = stmdif(:,si) * trlfac
            do ssi=1,Nstm
                stmdst(ssi) = sum(stmdif(:,ssi) * stmdif(:,ssi))
            end do
            minidx = minloc(stmdst)
            prdstm(ti,si) = minidx(1)

            conmtx(si, minidx(1)) = conmtx(si, minidx(1)) + prctrl
            
        end do
    end do

    deallocate(sumresp)
    deallocate(currsp)
    deallocate(stmdif)
    deallocate(stmdst)
    deallocate(minidx)

end subroutine mexFunction


subroutine bincount(x,c,n,m)
    implicit none

    integer, intent(in) :: n,m
    real(8), dimension(n), intent(in) :: x
    mwSize, dimension(0:m-1), intent(out) :: c

    integer :: i

    c = 0.0d0
    do i=1,n
        c(x(i)) = c(x(i)) + 1.0d0
    end do
end subroutine bincount
