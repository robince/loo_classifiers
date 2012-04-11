function [conmtx,info] = quadratic(data)
%LOO_CLASSIFIERS.QUADRATIC Quadratic discriminant analysis with LOOCV.
%   [CONMTX,INFO] = LOO_CLASSIFIERS.QUADRATIC(DATA) performs leave-one-out
%   cross validation using a quadratic discriminant classifier (full
%   covariance matrix estimated for each class).
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
Ntrl1 = Ntrl - 1;
Ntrl2 = Ntrl - 2;

prdstm = zeros(Ntrl,Ncls);
conmtx = zeros(Ncls,Ncls);

ftravg = mean(data, 3);
ftrsum = Ntrl*ftravg;

ftrcov = zeros(Nftr,Nftr,Ncls);
invftrcov = zeros(Nftr,Nftr,Ncls);
logdetcov = zeros(1,Ncls);
for ci=1:Ncls
    ftrcov(:,:,ci) = cov(squeeze(data(:,ci,:))');
    invftrcov(:,:,ci) = inv(ftrcov(:,:,ci));
    logdetcov(ci) = sum(log(svd(ftrcov(:,:,ci))));
end

prctrl = 100 / Ntrl;
invfac = (Ntrl - 2) / (Ntrl - 1);
detfac = -Ncls*log(invfac) - log(Ntrl);

for ci=1:Ncls 
    curftravg = ftravg;
    curlogdetcov = logdetcov;
    curinvftrcov = invftrcov;
    A = invftrcov(:,:,ci);
    for ti=1:Ntrl
        % define current vector
        curtrl = data(:,ci,ti);       
        % update mean
        curftravg(:,ci) = (ftrsum(:,ci) - curtrl) ./ Ntrl1;
        
        % subtract mean
        %demeantrlold = curtrl - ftravg(:,ci);
        demeantrlnew = curtrl - curftravg(:,ci);
        
        % flip sign since subtracting from original matrix
        % not adding as in SM formula
        u = -demeantrlnew;
        v = demeantrlnew;      
        
        % update covariance matrices using
        % Shannon-Morrison formula
        % http://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula
        num = A*u*v'*A;
        den = Ntrl + v'*A*u;
        curinvftrcov(:,:,ci) = invfac .* (A - num./den);
        
        % update determinant using matrix determinant lemma
        % http://en.wikipedia.org/wiki/Matrix_determinant_lemma
        curlogdetcov(ci) =  detfac + log(den) + logdetcov(ci);
        %curlogdetcov(ci) = sum(log(svd((Ntrl1.*ftrcov(:,:,ci) + u*v')./Ntrl2)));
        
        tmpdst = curtrl(:, ones(Ncls,1))-curftravg;
        stmdst = zeros(Ncls,1);
        for cci=1:Ncls
            stmdst(cci) = tmpdst(:,cci)'*curinvftrcov(:,:,cci)*tmpdst(:,cci) + curlogdetcov(cci);
        end
        
        [~, prdstm(ti,ci)] = min(stmdst); 
        conmtx(ci,prdstm(ti,ci)) = conmtx(ci,prdstm(ti,ci)) + prctrl;
    end % trials
end % classes

opts.method = 'dr';
opts.bias   = 'pt';
opts.btsp   = 0;
opts.nt     = Ntrl;
info = information(reshape(prdstm,[1 Ntrl Ncls]),opts,'I');

