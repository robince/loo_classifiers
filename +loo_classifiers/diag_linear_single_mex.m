function [conmtx,info] = diag_linear_single_mex(data)
%LOO_CLASSIFIERS.LINEAR Linear discriminant analysis with LOOCV.
%   [CONMTX,INFO] = LOO_CLASSIFIERS.LINEAR(DATA) performs leave-one-out
%   cross validation using a linear discriminant classifier (pooled
%   covariance matrix).
%
%   DATA should be a (Nftr, Ncls, Ntrl) 3D array. It is assumed that there 
%   are the same number of trials for each class, and class priors are 
%   uniform.
%
%   CONMTX is the (Ncls, Ncls) confusion matrix where the i,j entry gives
%   the percentage of class i trials which were classsified as class j.
%   mean(diag(CONMTX)) gives the average correct performance.
%
%   INFO is the information in the confusion matrix.
%
[Nftr, Ncls, Ntrl] = size(data);

% call mex version
[conmtx, prdstm] = loo_classifiers.diag_linear_single_core(data);
% 
% if nargin<2
%     opts.method = 'dr';
%     opts.bias   = 'pt';
%     opts.btsp   = 0;
%     opts.nt     = Ntrl;
%     info = information(reshape(prdstm,[1 Ntrl Ncls]),opts,'I');
% else
%     info = 0;
% end
info = 0;
