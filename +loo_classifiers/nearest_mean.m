function [conmtx,info] = nearest_mean(data)
% nearest_mean (features, classes, trials)

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

