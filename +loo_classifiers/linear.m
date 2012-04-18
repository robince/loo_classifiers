function [conmtx,info] = linear(data)
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
Ntottrl = Ntrl * Ncls;
Ntrl1 = Ntrl - 1;
Ntottrl1 = Ntottrl - 1;

prdstm = zeros(Ntrl,Ncls);
conmtx = zeros(Ncls,Ncls);

ftravg = mean(data, 3);
ftrsum = Ntrl*ftravg;
% pooled feature variance
demeandata = bsxfun(@minus, data, ftravg);

ftrcov = cov(demeandata(:,:)');
invftrcov = inv(ftrcov);
%ftrvar = (Ntottrl-1)*var(demeandata(:,:), [], 2);

prctrl = 100 / Ntrl;
invfac = (Ntottrl - 2) / (Ntottrl - 1);

% naming consistent with SM formula
A = invftrcov;
invfacA = invfac .* A;
for ci=1:Ncls 
    curftravg = ftravg;
    for ti=1:Ntrl
        % define current vector
        curtrl = data(:,ci,ti);       
        % update mean
        curftravg(:,ci) = (ftrsum(:,ci) - curtrl) ./ Ntrl1;
        % update covariance
        % subtract mean
        demeantrlold = curtrl - ftravg(:,ci);
        demeantrlnew = curtrl - curftravg(:,ci);
        
        % flip sign since subtracting from original matrix
        % not adding as in SM formula
        u = -demeantrlold;
        v = demeantrlnew;
        
        % update pooled covariance matrix using
        % Shannon-Morrison formula
        % http://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula
        num = (A*u)*(v'*A);
        den = (Ntottrl1 + u'*A*v) / invfac;
        %curinvftrcov = invfac .* (A - num./den);
        curinvftrcov = invfacA - (num./den);
        
        tmpdst = curtrl(:, ones(Ncls,1))-curftravg;       
%         stmdst = zeros(Ncls,1);
%         for cci=1:Ncls
%            stmdst(cci) = tmpdst(:,cci)'*curinvftrcov*tmpdst(:,cci); 
%         end
        % vectorized version
        stmdst = diag(tmpdst'*curinvftrcov*tmpdst);
        
        [~, prdstm(ti,ci)] = min(stmdst); 
        conmtx(ci,prdstm(ti,ci)) = conmtx(ci,prdstm(ti,ci)) + prctrl;
    end % trials
end % classes

opts.method = 'dr';
opts.bias   = 'pt';
opts.btsp   = 0;
opts.nt     = Ntrl;
info = information(reshape(prdstm,[1 Ntrl Ncls]),opts,'I');

