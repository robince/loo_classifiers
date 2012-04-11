function [conmtx,info] = diag_linear(data)
%LOO_CLASSIFIERS.DIAG_LINEAR Diagonal linear discriminant analysis with LOOCV.
%   [CONMTX,INFO] = LOO_CLASSIFIERS.DIAG_LINEAR(DATA) performs 
%   leave-one-out cross validation using a diagonal linear discriminant 
%   classifier (pooled diagonal covariance matrix).
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
Ntottrl = Ntrl * Ncls;

prdstm = zeros(Ntrl,Ncls);
conmtx = zeros(Ncls,Ncls);

ftravg = mean(data, 3);
ftrsum = Ntrl*ftravg;
% pooled feature variance
demeandata = bsxfun(@minus, data, ftravg);
ftrvar = (Ntottrl-1)*var(demeandata(:,:), [], 2);
prctrl = 100 / Ntrl;
% degrees of freedom for pooled variance with 1 trial out
dof = Ntottrl-1-Ncls;
Ntrl1 = Ntrl - 1;


for ci=1:Ncls 
    curftravg = ftravg;
    for ti=1:Ntrl
        % define current vector
        curtrl = data(:,ci,ti);       
        % update mean
        curftravg(:,ci) = (ftrsum(:,ci) - curtrl) ./ Ntrl1;
        % update variance using Knuth online algorithm
        % http://www.johndcook.com/standard_deviation.html
        % subtract mean
        demeantrlold = curtrl - ftravg(:,ci);
        demeantrlnew = curtrl - curftravg(:,ci);
        curftrvar = (ftrvar - demeantrlnew.*demeantrlold) / dof;
        
        %stmdst = sum( bsxfun(@rdivide, bsxfun(@minus, curtrl, curftravg).^2, curftrvar), 1);
        % much faster: (at least for small arrays)
        tmpdst = curtrl(:, ones(Ncls,1))-curftravg;
        stmdst = sum( (tmpdst.*tmpdst) ./ curftrvar(:, ones(Ncls,1)), 1); 
        
        [~, prdstm(ti,ci)] = min(stmdst); 
        conmtx(ci,prdstm(ti,ci)) = conmtx(ci,prdstm(ti,ci)) + prctrl;
    end % trials
end % classes

opts.method = 'dr';
opts.bias   = 'pt';
opts.btsp   = 0;
opts.nt     = Ntrl;
info = information(reshape(prdstm,[1 Ntrl Ncls]),opts,'I');

