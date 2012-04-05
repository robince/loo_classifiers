function [conmtx, info] = poisson_bayes(resp)
% Poisson Naive bayes decoding with leave-one-out cross validation
% DIY implementation
% ASSUMES INTEGER INPUT FROM 0
%
% Same syntax as Lin_decoding
% resp (datapoints, stimuli, trials)
[Nftr, Nstm, Ntrl] = size(resp);

% mean as max lik poison lamda
ftravg = mean(resp,3);
% handle zeros
% no spikes in training shouldnt imply resp with 1 spike there has 
% zero likelihood of being in that group.
zsmooth = 10;
if zsmooth, ftravg(ftravg==0) = 1 ./ (zsmooth*Ntrl); end

% cache factorial values for poisson likelihoods
maxcnt = max(resp(:));
facche = factorial(0:maxcnt)';
% cache exponentiation
expftravg = exp(-ftravg);


% confusion matrix
conmtx = zeros(Nstm,Nstm);
% predicted stimuli
prdstm = zeros(Ntrl,Nstm);
prctrl = 100 / Ntrl;

for si=1:Nstm
    curexpftravg = expftravg;
    curftravg = ftravg;
    for ti=1:Ntrl
        %curprb(:) = 0;
        thisr = resp(:,si,ti); 
        curftravg(:,si) = (Ntrl*ftravg(:,si) - thisr) ./ (Ntrl-1);
        curexpftravg(:,si) = exp(-curftravg(:,si));

        curlik = curexpftravg.*bsxfun(@rdivide,bsxfun(@power,curftravg,thisr),facche(thisr+1));

        curstmlik = prod(curlik);
        [~, prdstm(ti,si)] = max(curstmlik);

        conmtx(si, prdstm(ti,si)) = conmtx(si, prdstm(ti,si)) + prctrl;
    end
end


opts.method = 'dr';
opts.bias   = 'pt';
opts.btsp   = 0;
opts.nt     = Ntrl;
info = information(reshape(prdstm,[1 Ntrl Nstm]),opts,'I');
%info = info(1) - mean(info(2:end));

