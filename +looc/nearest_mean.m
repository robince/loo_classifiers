function [conmtx,I] = nearest_mean(X,Y,Ncls)
%LOO_CLASSIFIERS.NEAREST_MEAN Nearest mean / template matching with LOOCV.
%   [CONMTX,INFO] = LOO_CLASSIFIERS.NEAREST_MEAN(X,Y) performs leave-one-out
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
    error('nearest_mean: Class labels do not match data')
end
Y = Y(:);

prdY = zeros(Ntrl,1);
conmtx = zeros(Ncls,Ncls);

sum_data = zeros(Ncls, Nftr);
for fi=1:Nftr
    sum_data(:,fi) = accumarray(Y, X(:,fi));
end
Ntrlcls = accumarray(Y,1);

% switch to features first
data = X';
sum_data = sum_data';

prctrl = 100 ./ Ntrlcls;
trlratio = Ntrlcls ./ (Ntrlcls-1);
clsones = ones(Ncls, 1);
for ti=1:Ntrl
    thscls = Y(ti);
    % define current vector
    curtrl = data(:,ti) * Ntrlcls(thscls);
    tmpval = sum_data - curtrl(:, clsones);
    tmpval(:,thscls) = tmpval(:,thscls) * trlratio(thscls);
    clsdst = sum(tmpval .* tmpval, 1);
    [~, prdcls] = min(clsdst);
    conmtx(thscls,prdcls) = conmtx(thscls,prdcls) + prctrl(thscls);
    prdY(ti) = prdcls;
end
% [Y, prdY]
I = info.calc_info(Y-1,Ncls,prdY-1,Ncls,Ntrl);

