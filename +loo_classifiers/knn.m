function [conmtx,info] = knn(data,k,metric)
%LOO_CLASSIFIERS.KNN k-Nearest Neighbours with LOOCV.
%   [CONMTX,INFO] = LOO_CLASSIFIERS.KNN(DATA,K,METRIC) performs leave-one-out
%   cross validation using a k-nearest neighbours classifier.
%   Assigment is by majority rule with nearest point tie break.
%
%   DATA should be a (Nftr, Ncls, Ntrl) 3D array. It is assumed that there 
%   are the same number of trials for each class, and class priors are 
%   uniform.
%
%   K is the integer number of neighbours. METRIC is the distance metric 
%   to use.
%
%   CONMTX is the (Ncls, Ncls) confusion matrix where the i,j entry gives
%   the percentage of class i trials which were classsified as class j.
%   mean(diag(CONMTX)) gives the average correct performance.
%
%   INFO is the information in the confusion matrix.
%

[Nftr, Ncls, Ntrl] = size(data);

% build search object
nsmethod = 'kdtree';
%nsmethod = 'exhaustive';
NS = createns(data(:,:)','NSMethod',nsmethod,'Distance',metric);
cls = repmat(0:(Ncls-1),[1 Ntrl]);

prdstm = zeros(Ntrl,Ncls);
conmtx = zeros(Ncls,Ncls);

prctrl = 100 / Ntrl;

curidx = 1;
for ti=1:Ntrl
    for ci=1:Ncls
        % k+1 because we will obtain current point in results
        [idx, dst] = knnsearch(NS, data(:,ci,ti)','k',k+1,'IncludeTies',true);
        % remove current leave out
        % knncls indexes classes from 0 (for bincount)
        knncls = cls(idx{1}(idx{1}~=curidx));
        curidx = curidx + 1;
        % count votes
        votcnt = bincount(knncls,Ncls);
        [cnt, prestm] = max(votcnt);
        if cnt < k/2
            % might have a tie
            tieidx = find(votcnt==cnt);
            Ntie = length(tieidx);
            if Ntie>1
                % have a tie
                tiemin = zeros(1,length(tieidx));
                for tie=1:Ntie
                   tiemin(tie) = min(dst{1}(knncls == (tieidx(tie)-1)));
                end
                [~, nrst] = min(tiemin);
                prestm = tieidx(nrst);
            end % tie
        end % tie check
            
        prdstm(ti,ci) = prestm;
        conmtx(ci, prestm) = conmtx(ci, prestm) + prctrl;
    end % classes
end % trials

opts.method = 'dr';
opts.bias   = 'pt';
opts.btsp   = 0;
opts.nt     = Ntrl;
info = information(reshape(prdstm,[1 Ntrl Ncls]),opts,'I');

