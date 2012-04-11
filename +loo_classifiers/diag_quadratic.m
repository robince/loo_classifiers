function [conmtx,info] = diag_quadratic(data)
% diag_quadratic (features, classes, trials)

[Nftr, Ncls, Ntrl] = size(data);
Ntottrl = Ntrl * Ncls;
Ntrl1 = Ntrl - 1;
Ntrl2 = Ntrl - 2;

prdstm = zeros(Ntrl,Ncls);
conmtx = zeros(Ncls,Ncls);

ftravg = mean(data, 3);
ftrsum = Ntrl*ftravg;
ftrvar = var(data, [], 3);
prctrl = 100 / Ntrl;

for ci=1:Ncls 
    curftravg = ftravg;
    curftrvar = ftrvar;
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
        curftrvar(:,ci) = (Ntrl1.*ftrvar(:,ci) - demeantrlnew.*demeantrlold) / Ntrl2;
        curlogdet = sum(log(curftrvar),1);
        
        tmpdst = curtrl(:, ones(Ncls,1))-curftravg;
        stmdst = sum( (tmpdst.*tmpdst) ./ curftrvar, 1) + curlogdet; 
        
        [~, prdstm(ti,ci)] = min(stmdst); 
        conmtx(ci,prdstm(ti,ci)) = conmtx(ci,prdstm(ti,ci)) + prctrl;
    end % trials
end % classes

opts.method = 'dr';
opts.bias   = 'pt';
opts.btsp   = 0;
opts.nt     = Ntrl;
info = information(reshape(prdstm,[1 Ntrl Ncls]),opts,'I');

