% Need the Fortran 95 MatlabAPI Lite
% https://github.com/robince/MatlabAPI_Lite/


% mex options
DEBUG = false;
VERBOSE = false;
LARGEARRAY = true;

% path to matlab api
MAPI_INC = '/Users/robince/git/MatlabAPI_lite';
MAPI_LIB = MAPI_INC;

if ispc
    OBJEXT = 'obj';
else
    OBJEXT = 'o';
end

ARGS = {};
if LARGEARRAY,  ARGS{end+1} = '-largeArrayDims';
else            ARGS{end+1} = '-compatibleArrayDims'; end

if DEBUG,       ARGS{end+1} = '-g'; end
if VERBOSE,     ARGS{end+1} = '-v'; end
ARGS{end+1} = ['-I' MAPI_INC];

%%
% +looc_sorted
PKGDIR = '+looc_sorted';
PKGARGS = ARGS;
PKGARGS{end+1} = '-outdir'; PKGARGS{end+1} = PKGDIR;

% bincount
MEXARGS = PKGARGS;
MEXARGS{end+1} = fullfile(PKGDIR,'bincount.f');
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImx.' OBJEXT]);
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImex.' OBJEXT]);
mex(MEXARGS{:});

% nearest_mean
MEXARGS = PKGARGS;
MEXARGS{end+1} = fullfile(PKGDIR,'nearest_mean_core.f');
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImx.' OBJEXT]);
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImex.' OBJEXT]);
mex(MEXARGS{:});

% linear
MEXARGS = PKGARGS;
MEXARGS{end+1} = fullfile(PKGDIR,'linear_core.f')
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImx.' OBJEXT]);
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImex.' OBJEXT]);
mex(MEXARGS{:});


% diag_linear
MEXARGS = PKGARGS;
MEXARGS{end+1} = fullfile(PKGDIR,'diag_linear_core.f');
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImx.' OBJEXT]);
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImex.' OBJEXT]);
mex(MEXARGS{:});

% diag_linear_single
MEXARGS = PKGARGS;
MEXARGS{end+1} = fullfile(PKGDIR,'diag_linear_single_core.f');
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImx.' OBJEXT]);
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImex.' OBJEXT]);
mex(MEXARGS{:});

%%
%%
% +looc
PKGDIR = '+looc';
PKGARGS = ARGS;
PKGARGS{end+1} = '-outdir'; PKGARGS{end+1} = PKGDIR;

% nearest_mean
MEXARGS = PKGARGS;
MEXARGS{end+1} = fullfile(PKGDIR,'nearest_mean_core.f');
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImx.' OBJEXT]);
MEXARGS{end+1} = fullfile(MAPI_LIB,['MatlabAPImex.' OBJEXT]);
mex(MEXARGS{:});
