function [conmtx,info] = diag_linear(data)
% nearest_mean (features, classes, trials)

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
        % update variance
        demeantrl = curtrl - ftravg(:,ci);
        curftrvar = (ftrvar - demeantrl.*demeantrl) / dof;
        
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

