%
% MEX build script for windows
%
% Tested using intel compiler

% Need the Fortran 95 MatlabAPI
% http://www.mathworks.com/matlabcentral/fileexchange/25934
% https://github.com/robince/MatlabAPI/

% Path to MatlabAPI 
MAPI_INC = 'D:\robin\MatlabAPI';
MAPI_LIB = 'D:\robin\MatlabAPI';
PKGDIR = '+loo_classifiers';
LARGEARRAY = true;
DEBUG = false;
VERBOSE = false;

if ispc
    LIBEXT = '.obj';
else
    LIBEXT = '.lib';
end

ARGS = {};
if LARGEARRAY,  ARGS{end+1} = '-largeArrayDims';
else            ARGS{end+1} = '-compatibleArrayDims'; end

if DEBUG,       ARGS{end+1} = '-g'; end
if VERBOSE,     ARGS{end+1} = '-v'; end
ARGS{end+1} = '-outdir'; ARGS{end+1} = PKGDIR;
ARGS{end+1} = ['-I' MAPI_INC];
%%
% bincount
MEXARGS = ARGS;
MEXARGS{end+1} = fullfile(PKGDIR,'bincount.f');
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImx' LIBEXT]);
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImex' LIBEXT]);
mex(MEXARGS{:});

% nearest_mean
MEXARGS = ARGS;
MEXARGS{end+1} = fullfile(PKGDIR,'nearest_mean_core.f');
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImx' LIBEXT]);
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImex' LIBEXT]);
mex(MEXARGS{:});
%%
% linear
MEXARGS = ARGS;
% MEXARGS{end+1} = '-liomp5md'; % link omp dynamically as per docs
MEXARGS{end+1} = fullfile(PKGDIR,'linear_core.f');
%MEXARGS{end+1} = 'C:\Program Files (x86)\Intel\Composer XE 2011 SP1\mkl\lib\intel64\mkl_rt.lib';
MEXARGS{end+1} = 'mkl_lapack95_ilp64.lib';
MEXARGS{end+1} = 'mkl_blas95_ilp64.lib';
MEXARGS{end+1} = 'mkl_intel_ilp64.lib';
% MEXARGS{end+1} = 'mkl_intel_thread.lib';
MEXARGS{end+1} = 'mkl_sequential.lib';
MEXARGS{end+1} = 'mkl_core.lib';
% MEXARGS{end+1} = 'libiomp5md.lib';
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImx' LIBEXT]);
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImex' LIBEXT]);

mex(MEXARGS{:});
