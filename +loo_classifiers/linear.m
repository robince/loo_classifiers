function [conmtx,info] = linear(data)
% linear classifier (features, classes, trials)

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
        num = A*u*v'*A;
        dem = Ntottrl1 + u'*A*v;
        curinvftrcov = invfac .* (A - num./dem);
        
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

