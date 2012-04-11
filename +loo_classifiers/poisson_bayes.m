function [conmtx, info] = poisson_bayes(data)
%LOO_CLASSIFIERS.POISSON_BAYES Poisson Naive Bayes with LOOCV.
%   [CONMTX,INFO] = LOO_CLASSIFIERS.POISSON_BAYES(DATA) performs
%   leave-one-out cross validation using a Poisson naive bayes classifier
%   (class conditional feature probabilities are independent poisson).
%
%   DATA should be a (Nftr, Ncls, Ntrl) 3D array. It is assumed that there 
%   are the same number of trials for each class, and class priors are 
%   uniform. DATA MUST BE INTEGER COUNTS.
%
%   CONMTX is the (Ncls, Ncls) confusion matrix where the i,j entry gives
%   the percentage of class i trials which were classsified as class j.
%   mean(diag(CONMTX)) gives the average correct performance.
%
%   INFO is the information in the confusion matrix.
%

[Nftr, Ncls, Ntrl] = size(data);

% mean as max lik poison lamda
ftravg = mean(data,3);
% handle zeros
% no spikes in training shouldnt imply data with 1 spike there has 
% zero likelihood of being in that group.
zsmooth = 10;
if zsmooth, ftravg(ftravg==0) = 1 ./ (zsmooth*Ntrl); end

% cache factorial values for poisson likelihoods
maxcnt = max(data(:));
facche = factorial(0:maxcnt)';
% cache exponentiation
expftravg = exp(-ftravg);


% confusion matrix
conmtx = zeros(Ncls,Ncls);
% predicted stimuli
prdstm = zeros(Ntrl,Ncls);
prctrl = 100 / Ntrl;

for si=1:Ncls
    curexpftravg = expftravg;
    curftravg = ftravg;
    for ti=1:Ntrl
        %curprb(:) = 0;
        curtrl = data(:,si,ti);
        curftravg(:,si) = (Ntrl*ftravg(:,si) - curtrl) ./ (Ntrl-1);
        curexpftravg(:,si) = exp(-curftravg(:,si));

        %curlika = curexpftravg.*bsxfun(@rdivide,bsxfun(@power,curftravg,curtrl),facche(curtrl+1));
        % bsxfun is slow (for these array sizes)
        tmpnum = curftravg.^curtrl(:,ones(Ncls,1));
        tmpden = facche(curtrl+1);
        curlik = curexpftravg.*(tmpnum./tmpden(:,ones(Ncls,1)));
        
        curstmlik = prod(curlik);
        [~, prdstm(ti,si)] = max(curstmlik);

        conmtx(si, prdstm(ti,si)) = conmtx(si, prdstm(ti,si)) + prctrl;
    end
end


opts.method = 'dr';
opts.bias   = 'pt';
opts.btsp   = 0;
opts.nt     = Ntrl;
info = information(reshape(prdstm,[1 Ntrl Ncls]),opts,'I');
%info = info(1) - mean(info(2:end));

