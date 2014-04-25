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
    mwSize, parameter :: One=1
    mwPointer :: mxC

    if( nrhs /= 2 ) then
        call mexErrMsgTxt("This function takes 2 inputs")
    endif
    if( nlhs > 1) then
        call mexErrMsgTxt("This function returns 1 output")
    endif

    call fpGetPr( X, prhs(1) )
    call fpGetPr( M, prhs(2) )

    if( (.not.associated(X)) .or. (.not.associated(M)) ) then
        call mexErrMsgTxt("Problem with inputs: check types and dimensions")
    endif
    
    mi = M
    n = size(X)
    mxC = mxCreateNumericMatrix( mi, One, mxDOUBLE_CLASS, mxREAL )
    call fpGetPr( X, mxC )

    C = 0.0d0
    do i=1,n
        xi = x(i) + 1
        c(xi) = c(xi) + 1.0d0
    end do

    plhs(1) = mxC

end subroutine mexFunction


