function [conmtx,info] = multinom_bayes(data, lambda)
%LOO_CLASSIFIERS.MULTINOM_BAYES Multinomial Naive Bayes with LOOCV.
%   [CONMTX,INFO] = LOO_CLASSIFIERS.MULTINOM_BAYES(DATA,LAMBDA) performs
%   leave-one-out cross validation using a multinomial naive bayes classifier
%   (class conditional feature probabilities are independent multinomial).
%
%   DATA should be a (Nftr, Ncls, Ntrl) 3D array. It is assumed that there 
%   are the same number of trials for each class, and class priors are 
%   uniform. DATA MUST BE INTEGER COUNTS.
%
%   LAMBDA is the parameter for Laplace smoothing used when estimating
%   the probabilities.
%
%   CONMTX is the (Ncls, Ncls) confusion matrix where the i,j entry gives
%   the percentage of class i trials which were classsified as class j.
%   mean(diag(CONMTX)) gives the average correct performance.
%
%   INFO is the information in the confusion matrix.
%

[Nftr, Ncls, Ntrl] = size(data);

% dimension of each feature
ftrdim = max(data(:,:),[],2)+1;
% overall maximum feature dim
maxdim = max(ftrdim);

% counts for each feature value
ftrcnt = zeros(maxdim, Nftr, Ncls);

% build counts from full training data
for fi=1:Nftr
    for ci=1:Ncls
        ftrcnt(:,fi,ci) = loo_classifiers.bincount(data(fi,ci,:),maxdim);
    end
end
% smooth everywhere since non-occuring values will never
% be evaluated... only normalisation needs to keep track of 
% individual feature dimensions
ftrcnt = ftrcnt + lambda;

% normalise to probabilities
nrmval = Ntrl + lambda*ftrdim;
% normalise using dimension of each feature
ftrprb = bsxfun(@rdivide, ftrcnt, reshape(nrmval, [1 Nftr 1]));
% or normalise using largest dimension
%ftrprb = ftrcnt ./ (Ntrl + lambda*maxdim);

% build stim index for sub2ind 
stmidx = cell(Ncls,1);
for ci=1:Ncls
    stmidx{ci} = (ci-1)*ones(Nftr,1);
end
stmidx = cell2mat(stmidx);
yidx = repmat((0:(Nftr-1))', [Ncls 1]);
offsetidx = yidx*maxdim + stmidx*(maxdim*Nftr);

prdstm = zeros(Ntrl,Ncls);
conmtx = zeros(Ncls,Ncls);
prctrl = 100 / Ntrl;

%curprb = zeros(maxdim, Nftr);
for ci=1:Ncls
    curftrprb = ftrprb;
    for ti=1:Ntrl
        %curprb(:) = 0;
        % turns out for small sizes faster to allocate then reset
        curprb = zeros(maxdim, Nftr);
        curtrl = data(:,ci,ti)+1;
        for fi=1:Nftr
            curprb(curtrl(fi),fi) = 1 ./ (Ntrl + lambda*ftrdim(fi));
            %curprb(curtrl(fi),fi) = 1 ./ (Ntrl + lambda*maxdim);
        end
        % remove current trial from counts
        curftrprb(:,:,ci) = ftrprb(:,:,ci) - curprb;

        % normalise and predict
        nrmstm = nrmval ./ ((Ntrl-1) + lambda*ftrdim);
        nrmstm = nrmstm';
        curftrprb(:,:,ci) = curftrprb(:,:,ci) .* nrmstm(ones(maxdim,1),:);
        %curftrprb(:,:,ci) = curftrprb(:,:,ci) .* ((Ntrl+lambda*maxdim) ./ ((Ntrl-1)+lambda*maxdim)); 

        %xidx = repmat(curtrl, [Ncls 1]); 
        % repmat is slow
        xidx = curtrl(:,ones(1,Ncls));
        xidx = xidx(:);

        %idx = sub2ind(size(curftrprb), ...
                      %repmat(data(:,ci,ti)+1, [Ncls 1]), ...
                      %repmat((1:Nftr)', [Ncls 1]), ...
                      %stmidx);
        % sub2ind is slow
        idx = xidx + offsetidx;
        curlik = reshape(curftrprb(idx), [Nftr Ncls]);
        curstmlik = prod(curlik);

        [~, prdstm(ti,ci)] = max(curstmlik);
        conmtx(ci, prdstm(ti,ci)) = conmtx(ci, prdstm(ti,ci)) + prctrl;
    end % trials
end % classes


opts.method = 'dr';
opts.bias   = 'pt';
opts.btsp   = 0;
opts.nt     = Ntrl;
info = information(reshape(prdstm,[1 Ntrl Ncls]),opts,'I');

