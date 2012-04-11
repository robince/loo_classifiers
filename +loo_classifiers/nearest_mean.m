function [conmtx,info] = nearest_mean(data)
%LOO_CLASSIFIERS.NEAREST_MEAN Nearest mean / template matching with LOOCV.
%   [CONMTX,INFO] = LOO_CLASSIFIERS.NEAREST_MEAN(DATA) performs leave-one-out
%   cross validation using a nearest mean (template matching) discriminant 
%   classifier (pooled diagonal covariance with equal variances).
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

prdstm = zeros(Ntrl,Ncls);
conmtx = zeros(Ncls,Ncls);

sum_data = sum(data, 3);
prctrl = 100 / Ntrl;
data = data .* Ntrl;

for ci=1:Ncls 
    for ti=1:Ntrl
        % define current vector
        curtrl = data(:,ci,ti);       
        tmpval = sum_data - curtrl(:, ones(Ncls, 1));
        tmpval(:,ci) = tmpval(:,ci) * (Ntrl / (Ntrl-1));
        stmdst = sum(tmpval .* tmpval, 1);
        [~, prdstm(ti,ci)] = min(stmdst); 
        conmtx(ci,prdstm(ti,ci)) = conmtx(ci,prdstm(ti,ci)) + prctrl;
    end % trials
end % classes

opts.method = 'dr';
opts.bias   = 'pt';
opts.btsp   = 0;
opts.nt     = Ntrl;
info = information(reshape(prdstm,[1 Ntrl Ncls]),opts,'I');

