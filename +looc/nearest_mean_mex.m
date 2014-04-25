function [conmtx,I] = nearest_mean_mex(X,Y,Ncls)
%LOO_CLASSIFIERS.NEAREST_MEAN Nearest mean / template matching with LOOCV.
%   [CONMTX,INFO] = LOO_CLASSIFIERS.NEAREST_MEAN(X,Y,Ncls) performs leave-one-out
%   cross validation using a nearest mean (template matching) discriminant 
%   classifier (pooled diagonal covariance with equal variances).
%
%   X should be a (Ntrl, Nftr) matrix containing the feature data.
%   Y should be a length Ntrl vector labelling the class of each trial.
%   (integers indexed from 1)
%   Ncls is a scalar representing the number of classes.
%   Uniform class priors are used.
%
%   CONMTX is the (Ncls, Ncls) confusion matrix where the i,j entry gives
%   the percentage of class i trials which were classsified as class j.
%   mean(diag(CONMTX)) gives the average correct performance.
%
%   I is the information in the confusion matrix.
%
[Ntrl, Nftr] = size(X);
if length(Y) ~= Ntrl
    error('nearest_mean_mex: Class labels do not match data')
end

% call mex version
if ~strcmp(class(X), 'single')
    X = single(X);
end
if ~strcmp(class(Y), 'int16')
    Y = int16(Y);
end

[conmtx, prdY] = looc.nearest_mean_core(X',Y,Ncls);

I = info.calc_info(Y-1,Ncls,prdY-1,Ncls,Ntrl);
