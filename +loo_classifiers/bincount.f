!
! fast bincount
!
! DANGER - NO CHECKING AT ALL!
#include "fintrf.h"


subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use MatlabAPImex
    use MatlabAPImx
    implicit none
! ARG
    mwPointer plhs(*), prhs(*)
    integer(4) nlhs, nrhs
! LOC
    real(8), pointer :: X(:), M, C(:)
    mwSize :: i, mi, n, xi
    mwPointer :: mxC

    if( nrhs /= 2 ) then
        call mexErrMsgTxt("This function takes 2 inputs")
    endif
    if( nlhs > 1) then
        call mexErrMsgTxt("This function returns 1 output")
    endif

    X => fpGetPr1( prhs(1) );
    M => fpGetPr0( prhs(2) );

    if( (.not.associated(X)) .or. (.not.associated(M)) ) then
        call mexErrMsgTxt("Problem with inputs: check types and dimensions")
    endif
    
    mi = M
    n = size(X)
    C => fpAllocate1( mi )
    mxC = mxArrayHeader(C)

    C = 0.0d0
    do i=1,n
        xi = x(i) + 1
        c(xi) = c(xi) + 1.0d0
    end do

    plhs(1) = mxC

end subroutine mexFunction


